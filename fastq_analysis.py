"""
MinION Nanobody Analysis — v2
Ingest a full run folder — auto-detects all targets and runs the Test12.py
pipeline (trim+translate, easy-cluster, normalize) for each target group.
Results stored in SQLite. One ingest, all targets, instant reload.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
import shutil
import sqlite3
import subprocess
import tempfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components
from Bio import SeqIO
from Bio.Seq import Seq
from plotly.subplots import make_subplots


# ═══════════════════════════════════════════════════════════════════════════════
# ALL FUNCTIONS BELOW TAKEN DIRECTLY FROM Test12.py — DO NOT MODIFY LOGIC
# ═══════════════════════════════════════════════════════════════════════════════


# -----------------------------
# Helpers (folders / names)
# -----------------------------
def list_subfolders(data_dir: Path) -> list[str]:
    if not data_dir.exists():
        return []
    return sorted([p.name for p in data_dir.iterdir() if p.is_dir()])


def safe_project_name(name: str) -> str:
    name = (name or "").strip()
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "", name)
    safe = safe.strip(" .-_")
    return safe or "Project"


# -----------------------------
# Clipboard / export helpers
# -----------------------------
def clipboard_copy_button(
    text: str, label: str = "Copy nanobody sequences", key: str = "copy"
) -> None:
    """
    Renders an HTML button that copies text to the user's clipboard.
    Uses the browser Clipboard API; some environments may block clipboard access.
    """
    if text is None:
        text = ""

    # unique DOM ids per render to avoid collisions
    button_id = f"copybtn_{re.sub(r'[^A-Za-z0-9_]+', '', key)}{abs(hash(text)) % (10**9)}"
    msg_id = f"{button_id}_msg"
    text_json = json.dumps(text)

    components.html(
        f"""
        <div style="margin-top: 0.35rem; margin-bottom: 0.15rem;">
          <button id="{button_id}" style="
            background-color: rgb(14, 108, 255);
            color: white;
            border: none;
            padding: 0.45rem 0.75rem;
            border-radius: 0.25rem;
            cursor: pointer;
            font-size: 0.95rem;">
            {label}
          </button>
          <span id="{msg_id}" style="margin-left: 0.6rem; font-size: 0.95rem;"></span>
        </div>
        <script>
          const btn = document.getElementById("{button_id}");
          const msg = document.getElementById("{msg_id}");
          const text = {text_json};

          btn.addEventListener("click", async () => {{
            try {{
              await navigator.clipboard.writeText(text);
              msg.textContent = "Copied to clipboard.";
              setTimeout(() => msg.textContent = "", 2000);
            }} catch (e) {{
              msg.textContent = "Copy failed (browser blocked clipboard). Use the FASTA text below.";
              console.error(e);
            }}
          }});
        </script>
        """,
        height=60,
    )


def seq_df_to_fasta_text(seq_df: pd.DataFrame) -> str:
    """Convert a dataframe with columns ['cluster_head','aa_sequence'] into FASTA text."""
    if seq_df is None or seq_df.empty:
        return ""
    lines: list[str] = []
    for _, r in seq_df.iterrows():
        ch = str(r.get("cluster_head", "")).strip()
        seq = str(r.get("aa_sequence", "")).strip()
        if not ch or not seq:
            continue
        lines.append(f">{ch}\n{seq}")
    return ("\n".join(lines) + "\n") if lines else ""


def describe_file_for_log(p: Path) -> str:
    if p is None:
        return "(None)"
    try:
        p = Path(p)
    except Exception:
        return str(p)
    if not p.exists():
        return f"{p} [MISSING]"
    if p.is_dir():
        return f"{p} [DIR]"
    try:
        size = p.stat().st_size
        return f"{p} [OK, {size} bytes]"
    except Exception:
        return f"{p} [OK]"


# -----------------------------
# FASTQ helpers
# -----------------------------
def find_fastq_files(dir_path: Path) -> list[Path]:
    pats = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
    files: list[Path] = []
    for pat in pats:
        files.extend(sorted(dir_path.glob(f"*{pat}")))
    return files


def combine_fastqs(input_files: list[Path], combined_fastq: Path) -> None:
    """Concatenate FASTQ/FASTQ.GZ inputs into one plain FASTQ file."""
    combined_fastq.parent.mkdir(parents=True, exist_ok=True)
    with open(combined_fastq, "wb") as out_handle:
        for fp in input_files:
            if fp.suffix == ".gz":
                in_handle = gzip.open(fp, "rb")
            else:
                in_handle = open(fp, "rb")
            with in_handle:
                shutil.copyfileobj(in_handle, out_handle, length=1024 * 1024)  # 1 MB chunks


def iter_fastq_headers(fastq_path: Path):
    """
    Yields (read_header, read_id) for each FASTQ record.
    read_id is the first token before any whitespace.
    NOTE: Assumes 4-line FASTQ records (single-line seq/qual).
    """
    opener = gzip.open if fastq_path.suffix == ".gz" else open
    with opener(fastq_path, "rt", newline="") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            _seq = fh.readline()
            _plus = fh.readline()
            qual = fh.readline()
            if not qual:
                break

            header = header.rstrip("\n")
            header_no_at = header[1:] if header.startswith("@") else header
            read_id = header_no_at.split(" ", 1)[0]
            yield header_no_at, read_id


def count_fastq_records(files: list[Path]) -> int:
    """Count FASTQ records across a list of FASTQ/FASTQ.GZ files."""
    total = 0
    for fp in files:
        for _header, _read_id in iter_fastq_headers(fp):
            total += 1
    return total


# -----------------------------
# FASTQ -> protein FASTA
# -----------------------------
def trim_translate_fastq_to_fasta(
    input_fastq: Path,
    output_fasta: Path,
    START: str,
    END: str,
) -> tuple[int, int]:
    kept = 0
    discarded = 0
    START = (START or "").upper()
    END = (END or "").upper()

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(str(input_fastq), "fastq"):
            seq_str = str(record.seq).upper()

            start_idx = seq_str.find(START)
            end_idx = seq_str.find(END)

            if start_idx == -1 or end_idx == -1 or end_idx <= start_idx:
                discarded += 1
                continue

            sub_seq = seq_str[start_idx : end_idx + len(END)]
            codon_length = (len(sub_seq) // 3) * 3
            trimmed_seq = sub_seq[:codon_length]
            prot_seq = Seq(trimmed_seq).translate(to_stop=False)

            out_handle.write(f">{record.id}\n{prot_seq}\n")
            kept += 1

    return kept, discarded


def filter_aa_fasta(
    input_fasta: Path,
    output_fasta: Path,
    LENGTH: int,
) -> tuple[int, int]:
    kept = 0
    discarded = 0
    LENGTH = int(LENGTH)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(str(input_fasta), "fasta"):
            seq = str(record.seq)

            if len(seq) < LENGTH:
                discarded += 1
                continue
            if "*" in seq[:LENGTH]:
                discarded += 1
                continue

            SeqIO.write(record, out_handle, "fasta")
            kept += 1

    return kept, discarded


def fetch_fasta_sequences_by_id(fasta_path: Path, ids_in_order: list[str]) -> dict[str, str]:
    """Return {id: sequence} for the requested IDs (streaming parse; stops early when all found)."""
    wanted = set(map(str, ids_in_order))
    found: dict[str, str] = {}
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        rid = str(record.id)
        if rid in wanted and rid not in found:
            found[rid] = str(record.seq)
            if len(found) >= len(wanted):
                break
    return found


# -----------------------------
# MMseqs2 clustering
# -----------------------------
def run_mmseqs_easy_cluster(
    input_fasta: Path,
    output_prefix: Path,
    tmp_dir: Path,
    min_seq_id: float = 0.90,
    coverage: float = 0.9,
    cov_mode: int = 0,
    use_linclust: bool = False,
) -> tuple[list[str], str, str, list[Path]]:
    """
    Runs:
      mmseqs easy-cluster input.fasta out_prefix tmp_dir --min-seq-id 0.90 -c 0.9 --cov-mode 0
    Returns: (cmd, stdout, stderr, output_files)
    """
    if shutil.which("mmseqs") is None:
        raise RuntimeError(
            "mmseqs was not found in your PATH. Install MMseqs2 and ensure the mmseqs executable is available."
        )
    if not input_fasta.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    cluster_mode = "easy-linclust" if use_linclust else "easy-cluster"
    cmd = [
        "mmseqs",
        cluster_mode,
        str(input_fasta),
        str(output_prefix),
        str(tmp_dir),
        "--min-seq-id",
        str(min_seq_id),
        "-c",
        str(coverage),
        "--cov-mode",
        str(cov_mode),
    ]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "MMseqs2 failed.\n"
            f"Command: {' '.join(cmd)}\n\n"
            f"STDOUT:\n{proc.stdout}\n\n"
            f"STDERR:\n{proc.stderr}"
        )

    expected = [
        output_prefix.parent / f"{output_prefix.name}_cluster.tsv",
        output_prefix.parent / f"{output_prefix.name}_rep_seq.fasta",
        output_prefix.parent / f"{output_prefix.name}_all_seqs.fasta",
    ]
    out_files = [p for p in expected if p.exists()]
    if not out_files:
        out_files = sorted(output_prefix.parent.glob(output_prefix.name + "*"))

    return cmd, proc.stdout, proc.stderr, out_files


# -----------------------------
# cluster abundance + differential abundance
# -----------------------------
def build_cluster_counts_and_matrix(
    library_paths: dict[str, Path],
    cluster_tsv: Path,
    drop_unclustered: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Reads mmseqs cluster TSV (cluster_head, read_id), then scans original FASTQs per library
    and counts reads per (library_id, cluster_head).

    If drop_unclustered=False, reads not in the TSV are counted under cluster_head="__unclustered__".
    """
    if not cluster_tsv.exists():
        raise FileNotFoundError(f"Cluster TSV not found: {cluster_tsv}")

    cluster_all = pd.read_csv(
        cluster_tsv,
        sep="\t",
        header=None,
        names=["cluster_head", "read_id"],
        dtype=str,
    )

    # Map member read_id -> cluster_head
    read_to_cluster = dict(
        zip(cluster_all["read_id"].astype(str), cluster_all["cluster_head"].astype(str))
    )

    counts = Counter()

    for lib_id, lib_path in library_paths.items():
        lib_path = Path(lib_path)
        for fq in find_fastq_files(lib_path):
            for _read_header, read_id in iter_fastq_headers(fq):
                ch = read_to_cluster.get(read_id)
                if ch is None:
                    if not drop_unclustered:
                        counts[(lib_id, "__unclustered__")] += 1
                    continue
                counts[(lib_id, ch)] += 1

    cluster_counts = pd.DataFrame(
        [{"library_id": lib, "cluster_head": ch, "n": n} for (lib, ch), n in counts.items()]
    )

    if cluster_counts.empty:
        count_matrix = pd.DataFrame({"cluster_head": []})
        return cluster_counts, count_matrix

    count_matrix = (
        cluster_counts.pivot_table(
            index="cluster_head",
            columns="library_id",
            values="n",
            fill_value=0,
            aggfunc="sum",
        )
        .reset_index()
    )

    return cluster_counts, count_matrix


def plot_library_abundance_distribution(
    cluster_counts: pd.DataFrame,
    out_pdf: Path,
    out_tsv: Path,
) -> tuple[plt.Figure, pd.DataFrame]:
    plot_rows = []
    if not cluster_counts.empty:
        cluster_counts = cluster_counts.copy()
        cluster_counts = cluster_counts[cluster_counts["n"] > 0]
        cluster_counts["log10_counts"] = np.log10(cluster_counts["n"].astype(float))

    for lib_id, df_lib in cluster_counts.groupby("library_id"):
        x = df_lib["log10_counts"].to_numpy()
        if x.size == 0:
            continue

        counts, edges = np.histogram(x, bins="sturges")
        mids = (edges[:-1] + edges[1:]) / 2

        y = np.full_like(mids, np.nan, dtype=float)
        mask = counts > 0
        y[mask] = np.log10(counts[mask].astype(float))

        plot_rows.append(pd.DataFrame({"library_id": lib_id, "x": mids, "y": y}))

    plot_data_abundance = (
        pd.concat(plot_rows, ignore_index=True)
        if plot_rows
        else pd.DataFrame({"library_id": [], "x": [], "y": []})
    )

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except Exception:
        pass

    fig, ax = plt.subplots(figsize=(5, 5))
    if not plot_data_abundance.empty:
        for lib_id, df_lib in plot_data_abundance.groupby("library_id"):
            ax.plot(df_lib["x"], df_lib["y"], linewidth=1.9, label=lib_id)

    ax.set_title("Abundance Distribution by Library")
    ax.set_xlabel("Abundance")
    ax.set_ylabel("N Nanobody Sequence Clusters")
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xticklabels(["0", "10", "100", "1000", "10000"])
    ax.set_yticks([0, 1, 2, 3, 4, 5])
    ax.set_yticklabels(["0", "1e1", "1e2", "1e3", "1e4", "1e5"])
    ax.legend(title="Library ID", loc="upper left", bbox_to_anchor=(0.62, 0.86))
    fig.tight_layout()

    fig.savefig(out_pdf)
    plot_data_abundance.to_csv(out_tsv, sep="\t", index=False)

    return fig, plot_data_abundance


def plot_abundance_vs_differential(
    count_matrix: pd.DataFrame,
    baseline_lib: str,
    condition_libs: list[str],
    baseline_threshold: int,
    enrichment_threshold: float,
    out_pdf: Path,
    out_tsv: Path,
) -> tuple[plt.Figure, pd.DataFrame]:
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Ensure expected columns exist
    for col in [baseline_lib] + condition_libs:
        if col not in count_matrix.columns:
            count_matrix[col] = 0

    frames = []
    for cond2 in condition_libs:
        df = count_matrix.assign(
            condition1=baseline_lib,
            condition2=cond2,
            count1=count_matrix[baseline_lib].astype(float),
            count2=count_matrix[cond2].astype(float),
        )
        df["ratio"] = df["count2"] / df["count1"].replace(0, np.nan)
        frames.append(df)

    plot_data_diff = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if plot_data_diff.empty:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.set_title("Differential Abundance Analysis (no data)")
        fig.tight_layout()
        fig.savefig(out_pdf)
        plot_data_diff.to_csv(out_tsv, sep="\t", index=False)
        return fig, plot_data_diff

    plot_data_diff = plot_data_diff[
        np.isfinite(plot_data_diff["ratio"]) & (plot_data_diff["ratio"] > 0)
    ].copy()
    plot_data_diff["highlight"] = (
        (plot_data_diff["count1"] >= float(baseline_threshold))
        & (plot_data_diff["ratio"] >= float(enrichment_threshold))
    )

    # Significant-in-first-condition gets blue across both facets
    sig_ref_condition = condition_libs[0] if condition_libs else None
    sig_in_ref = (
        set(
            plot_data_diff.loc[
                (plot_data_diff["condition2"] == sig_ref_condition)
                & (plot_data_diff["highlight"]),
                "cluster_head",
            ]
        )
        if sig_ref_condition
        else set()
    )

    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except Exception:
        pass

    n_facets = max(1, len(condition_libs))
    fig, axes = plt.subplots(1, n_facets, figsize=(7, 5), sharex=True, sharey=True)
    if n_facets == 1:
        axes = [axes]

    for ax, cond2 in zip(axes, condition_libs):
        d = plot_data_diff[plot_data_diff["condition2"] == cond2]

        ax.axvline(np.log10(float(baseline_threshold)), color="gray")
        ax.axhline(np.log10(float(enrichment_threshold)), color="gray")

        d0 = d[~d["highlight"]]
        ax.scatter(
            np.log10(d0["count1"]),
            np.log10(d0["ratio"]),
            c="black",
            alpha=0.6,
            s=6,
            marker="o",
            linewidths=0,
            zorder=1,
        )

        d_blue = d[d["cluster_head"].isin(sig_in_ref)]
        ax.scatter(
            np.log10(d_blue["count1"]),
            np.log10(d_blue["ratio"]),
            c="blue",
            alpha=0.9,
            s=22,
            marker="o",
            linewidths=0,
            zorder=2,
        )

        # Red only for the *second* facet
        if len(condition_libs) >= 2 and cond2 == condition_libs[1]:
            d_red = d[d["highlight"]]
            ax.scatter(
                np.log10(d_red["count1"]),
                np.log10(d_red["ratio"]),
                c="darkred",
                alpha=0.9,
                s=22,
                marker="o",
                linewidths=0,
                zorder=3,
            )

        ax.set_title(cond2)
        ax.set_xlabel("Baseline Abundance")

    axes[0].set_ylabel("(Condition / Baseline) Abundance")

    for ax in axes:
        ax.set_xticks(np.log10([1, 10, 100, 1000]))
        ax.set_xticklabels(["1", "10", "100", "1000"])
        ax.set_yticks(np.log10([0.01, 1, 100, 10000]))
        ax.set_yticklabels(["1:100", "1:1", "100:1", "10000:1"])

    fig.suptitle(
        "Differential Abundance Analysis\n"
        f"Significance: Baseline Abundance >= {baseline_threshold}, Enrichment >= {enrichment_threshold} fold",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(out_pdf)
    plot_data_diff.to_csv(out_tsv, sep="\t", index=False)

    return fig, plot_data_diff


# NEW: Plotly version for interactive hover (ID + AA sequence)
def plot_abundance_vs_differential_plotly(
    count_matrix: pd.DataFrame,
    baseline_lib: str,
    condition_libs: list[str],
    baseline_threshold: int,
    enrichment_threshold: float,
    aa_seq_map: dict[str, str] | None = None,  # {cluster_head: aa_sequence}
) -> go.Figure:
    # Ensure expected columns exist
    for col in [baseline_lib] + condition_libs:
        if col not in count_matrix.columns:
            count_matrix[col] = 0

    frames = []
    for cond2 in condition_libs:
        df = count_matrix.assign(
            condition1=baseline_lib,
            condition2=cond2,
            count1=count_matrix[baseline_lib].astype(float),
            count2=count_matrix[cond2].astype(float),
        )
        df["ratio"] = df["count2"] / df["count1"].replace(0, np.nan)
        frames.append(df)

    plot_data_diff = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if plot_data_diff.empty:
        return go.Figure()

    plot_data_diff = plot_data_diff[
        np.isfinite(plot_data_diff["ratio"]) & (plot_data_diff["ratio"] > 0)
    ].copy()
    plot_data_diff["highlight"] = (
        (plot_data_diff["count1"] >= float(baseline_threshold))
        & (plot_data_diff["ratio"] >= float(enrichment_threshold))
    )

    # Significant-in-first-condition gets blue across both facets
    sig_ref_condition = condition_libs[0] if condition_libs else None
    sig_in_ref = (
        set(
            plot_data_diff.loc[
                (plot_data_diff["condition2"] == sig_ref_condition)
                & (plot_data_diff["highlight"]),
                "cluster_head",
            ]
        )
        if sig_ref_condition
        else set()
    )

    if aa_seq_map is None:
        aa_seq_map = {}
    plot_data_diff["aa_sequence"] = (
        plot_data_diff["cluster_head"].astype(str).map(aa_seq_map).fillna("")
    )

    # log coords
    plot_data_diff["x"] = np.log10(plot_data_diff["count1"].astype(float))
    plot_data_diff["y"] = np.log10(plot_data_diff["ratio"].astype(float))

    n_facets = max(1, len(condition_libs))
    fig = make_subplots(
        rows=1,
        cols=n_facets,
        subplot_titles=condition_libs if condition_libs else [""],
        shared_xaxes=True,
        shared_yaxes=True,
        horizontal_spacing=0.08,
    )

    x_vline = np.log10(float(baseline_threshold))
    y_hline = np.log10(float(enrichment_threshold))

    # Hover: show ID + AA sequence (plus numeric context)
    hovertemplate = (
        "Nanobody ID: %{customdata[0]}<br>"
        "AA sequence: %{customdata[1]}<br>"
        "Baseline count: %{customdata[2]:.0f}<br>"
        "Fold enrichment: %{customdata[3]:.3g}×"
        "<extra></extra>"
    )

    for i, cond2 in enumerate(condition_libs, start=1):
        d = plot_data_diff[plot_data_diff["condition2"] == cond2].copy()

        d0 = d[~d["highlight"]]
        d_blue = d[d["cluster_head"].isin(sig_in_ref)]
        d_red = (
            d[d["highlight"]]
            if (len(condition_libs) >= 2 and cond2 == condition_libs[1])
            else d.iloc[0:0]
        )

        def add_scatter(dd: pd.DataFrame, color: str, size: int, name: str, showlegend: bool):
            if dd.empty:
                return
            customdata = np.stack(
                [
                    dd["cluster_head"].astype(str).to_numpy(),
                    dd["aa_sequence"].astype(str).to_numpy(),
                    dd["count1"].to_numpy(),
                    dd["ratio"].to_numpy(),
                ],
                axis=1,
            )
            fig.add_trace(
                go.Scatter(
                    x=dd["x"],
                    y=dd["y"],
                    mode="markers",
                    marker=dict(color=color, size=size, opacity=0.85),
                    name=name,
                    showlegend=showlegend,
                    customdata=customdata,
                    hovertemplate=hovertemplate,
                ),
                row=1,
                col=i,
            )

        add_scatter(d0, "black", 5, "not significant", showlegend=(i == 1))
        add_scatter(d_blue, "blue", 9, "sig in 1xpanned (blue)", showlegend=(i == 1))
        add_scatter(d_red, "darkred", 9, "sig in 2xpanned (red)", showlegend=(i == 1))

        fig.add_vline(x=x_vline, line_color="gray", line_width=1, row=1, col=i)
        fig.add_hline(y=y_hline, line_color="gray", line_width=1, row=1, col=i)

        fig.update_xaxes(
            title_text="Baseline Abundance",
            tickmode="array",
            tickvals=np.log10([1, 10, 100, 1000]),
            ticktext=["1", "10", "100", "1000"],
            row=1,
            col=i,
        )

    fig.update_yaxes(
        title_text="(Condition / Baseline) Abundance",
        tickmode="array",
        tickvals=np.log10([0.01, 1, 100, 10000]),
        ticktext=["1:100", "1:1", "100:1", "10000:1"],
        row=1,
        col=1,
    )

    fig.update_layout(
        title=(
            "Differential Abundance Analysis<br>"
            f"<sup>Significance: Baseline ≥ {baseline_threshold}, Enrichment ≥ {enrichment_threshold}×</sup>"
        ),
        height=520,
        margin=dict(l=60, r=20, t=80, b=60),
    )
    return fig


# -----------------------------
# Streamlit UI
# -----------------------------


# ═══════════════════════════════════════════════════════════════════════════════
# BARCODE FOLDER PARSING
# ═══════════════════════════════════════════════════════════════════════════════

def parse_barcode_label(folder_name: str) -> dict:
    name = folder_name.strip()
    bc_match = re.match(r"^(barcode\d+)", name, re.IGNORECASE)
    if not bc_match:
        return {}
    barcode = bc_match.group(1).lower()
    rest = name[len(barcode):].lstrip("_")
    if not rest:
        return {"barcode": barcode, "target": "unknown", "round": 1, "label": name}
    if re.search(r"_TG1$", rest, re.IGNORECASE) or re.match(r"^TG1$", rest, re.IGNORECASE):
        target = re.sub(r"_?TG1$", "", rest, flags=re.IGNORECASE).strip("_") or "TG1"
        return {"barcode": barcode, "target": target, "round": 0, "label": name}
    round_match = re.search(r"_R(\d+)$", rest, re.IGNORECASE)
    if round_match:
        round_num = int(round_match.group(1))
        target = rest[: round_match.start()].strip("_")
        return {"barcode": barcode, "target": target, "round": round_num, "label": name}
    return {"barcode": barcode, "target": rest, "round": 1, "label": name}


def discover_barcodes(folder: Path) -> List[dict]:
    """Discover all barcode folders in a run directory."""
    results = []
    scan_dir = folder / "fastq_pass" if (folder / "fastq_pass").exists() else folder
    for d in sorted(scan_dir.iterdir()):
        if not d.is_dir() or "unclassified" in d.name.lower():
            continue
        meta = parse_barcode_label(d.name)
        if not meta:
            continue
        files = find_fastq_files(d)
        if not files:
            continue
        meta["files"] = files
        meta["folder_path"] = str(d.resolve())
        results.append(meta)
    return results


def group_barcodes_by_target(barcodes: List[dict]) -> Dict[str, dict]:
    """
    Group barcodes by target. For each target return:
      {target: {tg1: bc_info, rounds: {1: bc_info, 2: bc_info, ...}}}
    TG1 (round=0) matching uses prefix match for cases like C-auris_Nm -> C-auris_Nm_CW.
    """
    groups: Dict[str, dict] = {}

    # First pass: collect all targets and their barcodes
    for bc in barcodes:
        target = bc.get("target", "unknown")
        round_num = bc.get("round", 1)
        if target not in groups:
            groups[target] = {"tg1": None, "rounds": {}}
        if round_num == 0:
            groups[target]["tg1"] = bc
        else:
            groups[target]["rounds"][round_num] = bc

    # Second pass: for targets with no TG1, try prefix match
    tg1_barcodes = [bc for bc in barcodes if bc.get("round") == 0]
    for target, grp in groups.items():
        if grp["tg1"] is None:
            for tg1_bc in tg1_barcodes:
                tg1_target = tg1_bc.get("target", "")
                if target.startswith(tg1_target):
                    grp["tg1"] = tg1_bc
                    break

    # Remove targets with no rounds (pure TG1 barcodes with no panning)
    groups = {t: g for t, g in groups.items() if g["rounds"]}

    return groups


# ═══════════════════════════════════════════════════════════════════════════════
# DATABASE
# ═══════════════════════════════════════════════════════════════════════════════

DEFAULT_DB                   = "~/minion_nanobody_v2.db"
DEFAULT_START                = "ATGGCC"
DEFAULT_END                  = "GGCGCGC"
DEFAULT_LENGTH               = 80
DEFAULT_BASELINE_THRESHOLD   = 100
DEFAULT_ENRICHMENT_THRESHOLD = 10.0

SCHEMA_SQL = """
PRAGMA journal_mode=WAL;

CREATE TABLE IF NOT EXISTS run (
  run_id          TEXT PRIMARY KEY,
  flowcell_id     TEXT,
  sample_id       TEXT,
  folder_path     TEXT,
  start_anchor    TEXT,
  end_anchor      TEXT,
  min_aa_len      INTEGER,
  mm_min_seq_id   REAL,
  mm_coverage     REAL,
  mm_cov_mode     INTEGER,
  use_linclust    INTEGER,
  n_targets       INTEGER,
  ingested_at     TEXT
);

CREATE TABLE IF NOT EXISTS target_run (
  run_id          TEXT,
  target          TEXT,
  control_path    TEXT,
  onex_path       TEXT,
  twox_path       TEXT,
  total_reads_control    INTEGER,
  total_reads_1xpanned   INTEGER,
  total_reads_2xpanned   INTEGER,
  kept_trimmed    INTEGER,
  kept_filtered   INTEGER,
  n_clusters      INTEGER,
  PRIMARY KEY (run_id, target)
);

CREATE TABLE IF NOT EXISTS cluster_counts (
  run_id        TEXT,
  target        TEXT,
  cluster_head  TEXT,
  library_id    TEXT,
  raw_count     INTEGER,
  norm_count    REAL,
  PRIMARY KEY (run_id, target, cluster_head, library_id)
);

CREATE TABLE IF NOT EXISTS cluster_sequences (
  run_id        TEXT,
  target        TEXT,
  cluster_head  TEXT,
  aa_sequence   TEXT,
  aa_length     INTEGER,
  PRIMARY KEY (run_id, target, cluster_head)
);
"""


def connect_db(db_path: Path) -> sqlite3.Connection:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path), timeout=30, check_same_thread=False)
    conn.execute("PRAGMA foreign_keys=ON;")
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA busy_timeout=10000;")
    conn.executescript(SCHEMA_SQL)
    return conn


def sql_df(conn: sqlite3.Connection, sql: str, params=()) -> pd.DataFrame:
    try:
        return pd.read_sql_query(sql, conn, params=params)
    except Exception:
        return pd.DataFrame()


# ═══════════════════════════════════════════════════════════════════════════════
# INGEST — whole run, one target at a time using Test12.py pipeline
# ═══════════════════════════════════════════════════════════════════════════════

def ingest_target(
    conn: sqlite3.Connection,
    run_id: str,
    target: str,
    control_dir: Path,
    onex_dir: Optional[Path],
    twox_dir: Optional[Path],
    START: str,
    END: str,
    LENGTH: int,
    mm_min_seq_id: float,
    mm_coverage: float,
    mm_cov_mode: int,
    use_linclust: bool,
    drop_unclustered: bool,
    progress_cb=None,
):
    """Run Test12.py pipeline for one target group and store results."""
    def log(msg):
        if progress_cb:
            progress_cb(f"  [{target}] {msg}")

    control_files = find_fastq_files(control_dir) if control_dir else []
    one_x_files   = find_fastq_files(onex_dir)    if onex_dir   else []
    two_x_files   = find_fastq_files(twox_dir)    if twox_dir   else []

    # If no TG1, use R1 as control and R2 as 1xpanned
    if not control_files and one_x_files and two_x_files:
        log("No TG1 — using R1 as control, R2 as 1xpanned")
        control_dir   = onex_dir
        control_files = one_x_files
        onex_dir      = twox_dir
        one_x_files   = two_x_files
        twox_dir      = None
        two_x_files   = []
    elif not control_files and one_x_files:
        log("No TG1 and only R1 — skipping (need at least 2 rounds)")
        return

    all_inputs = control_files + one_x_files + two_x_files
    if not all_inputs:
        log("No FASTQ files found — skipping")
        return

    with tempfile.TemporaryDirectory(prefix=f"nb_{target}_") as tmp:
        tmp_path       = Path(tmp)
        combined_fastq = tmp_path / "combined.fastq"
        trimmed_fasta  = tmp_path / "trimmed_proteins.fasta"
        filtered_fasta = tmp_path / "filtered.fasta"
        cluster_prefix = tmp_path / "clusters"
        cluster_tmp    = tmp_path / "mmseqs_tmp"
        cluster_tsv    = tmp_path / "clusters_cluster.tsv"

        # Count reads (Test12.py)
        log("Counting reads...")
        total_reads_control  = count_fastq_records(control_files)
        total_reads_1xpanned = count_fastq_records(one_x_files)
        total_reads_2xpanned = count_fastq_records(two_x_files)
        log(f"control={total_reads_control:,}  1x={total_reads_1xpanned:,}  2x={total_reads_2xpanned:,}")

        # Combine (Test12.py)
        log("Combining FASTQs...")
        combine_fastqs(all_inputs, combined_fastq)

        # Trim + translate (Test12.py)
        log(f"Trimming (START={START}, END={END}) and translating...")
        kept5, discarded5 = trim_translate_fastq_to_fasta(combined_fastq, trimmed_fasta, START=START, END=END)
        log(f"Kept {kept5:,}, discarded {discarded5:,}")

        # Filter (Test12.py)
        log(f"Filtering (min AA length={LENGTH})...")
        kept6, discarded6 = filter_aa_fasta(trimmed_fasta, filtered_fasta, LENGTH=int(LENGTH))
        log(f"Kept {kept6:,}, discarded {discarded6:,}")

        if kept6 <= 0:
            log("No sequences passed filtering — skipping. Check START/END anchors.")
            return

        # Cluster (Test12.py)
        mode_str = "easy-linclust" if use_linclust else "easy-cluster"
        log(f"Clustering with MMseqs2 {mode_str}...")
        run_mmseqs_easy_cluster(
            filtered_fasta, cluster_prefix, cluster_tmp,
            min_seq_id=mm_min_seq_id, coverage=mm_coverage,
            cov_mode=mm_cov_mode, use_linclust=use_linclust,
        )
        if not cluster_tsv.exists():
            log("Cluster TSV not found — skipping")
            return
        log("Clustering complete")

        # Count matrix (Test12.py)
        log("Building count matrix...")
        lib_paths = {"control": control_dir}
        if onex_dir:  lib_paths["1xpanned"] = onex_dir
        if twox_dir:  lib_paths["2xpanned"] = twox_dir
        cluster_counts_df, count_matrix = build_cluster_counts_and_matrix(
            library_paths=lib_paths,
            cluster_tsv=cluster_tsv,
            drop_unclustered=drop_unclustered,
        )
        for col in ["control", "1xpanned", "2xpanned"]:
            if col not in count_matrix.columns:
                count_matrix[col] = 0

        # Normalize (Test12.py — exact formula)
        log("Normalizing...")
        count_matrix["control_norm"] = count_matrix["control"].astype(float).round(0).astype("Int64")
        if total_reads_control > 0 and total_reads_1xpanned > 0:
            count_matrix["1xpanned_norm"] = (
                count_matrix["1xpanned"].astype(float)
                * (float(total_reads_control) / float(total_reads_1xpanned))
            ).round(0).astype("Int64")
        else:
            count_matrix["1xpanned_norm"] = pd.Series([pd.NA]*len(count_matrix), dtype="Int64")
        if total_reads_control > 0 and total_reads_2xpanned > 0:
            count_matrix["2xpanned_norm"] = (
                count_matrix["2xpanned"].astype(float)
                * (float(total_reads_control) / float(total_reads_2xpanned))
            ).round(0).astype("Int64")
        else:
            count_matrix["2xpanned_norm"] = pd.Series([pd.NA]*len(count_matrix), dtype="Int64")

        count_matrix = count_matrix.sort_values(
            by="2xpanned_norm", ascending=False, kind="mergesort", na_position="last"
        ).reset_index(drop=True)

        n_clusters = len(count_matrix)
        log(f"{n_clusters:,} clusters")

        # Fetch sequences (Test12.py)
        log("Fetching AA sequences...")
        all_heads = count_matrix["cluster_head"].astype(str).tolist()
        seq_map = fetch_fasta_sequences_by_id(filtered_fasta, all_heads)

        # Store
        log("Saving to database...")
        with conn:
            conn.execute(
                """INSERT OR REPLACE INTO target_run
                   (run_id, target, control_path, onex_path, twox_path,
                    total_reads_control, total_reads_1xpanned, total_reads_2xpanned,
                    kept_trimmed, kept_filtered, n_clusters)
                   VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
                (run_id, target,
                 str(control_dir) if control_dir else None,
                 str(onex_dir)    if onex_dir    else None,
                 str(twox_dir)    if twox_dir    else None,
                 total_reads_control, total_reads_1xpanned, total_reads_2xpanned,
                 kept5, kept6, n_clusters)
            )
            count_rows = []
            for _, row in count_matrix.iterrows():
                ch = str(row["cluster_head"])
                for lib, norm_col in [("control","control_norm"),("1xpanned","1xpanned_norm"),("2xpanned","2xpanned_norm")]:
                    raw = int(row.get(lib, 0) or 0)
                    nv  = row.get(norm_col, raw)
                    norm = float(nv) if nv is not pd.NA and nv is not None else float(raw)
                    count_rows.append((run_id, target, ch, lib, raw, norm))
            conn.executemany(
                "INSERT OR REPLACE INTO cluster_counts (run_id,target,cluster_head,library_id,raw_count,norm_count) VALUES (?,?,?,?,?,?)",
                count_rows
            )
            seq_rows = [(run_id, target, ch, seq, len(seq)) for ch, seq in seq_map.items()]
            conn.executemany(
                "INSERT OR IGNORE INTO cluster_sequences (run_id,target,cluster_head,aa_sequence,aa_length) VALUES (?,?,?,?,?)",
                seq_rows
            )
        log(f"Done ✅")


def ingest_run(
    fastq_dir: Path,
    db_path: Path,
    START: str,
    END: str,
    LENGTH: int,
    mm_min_seq_id: float,
    mm_coverage: float,
    mm_cov_mode: int,
    use_linclust: bool,
    drop_unclustered: bool,
    ext_tg1_dir: Optional[Path] = None,
    progress_cb=None,
) -> str:
    def log(msg):
        if progress_cb:
            progress_cb(msg)

    log("Scanning barcode folders...")
    barcodes = discover_barcodes(fastq_dir)
    if not barcodes:
        raise ValueError("No barcode folders found.")

    # Get run metadata from first FASTQ header
    run_id = flowcell_id = sample_id = "unknown"
    first_file = barcodes[0]["files"][0]
    opener = gzip.open if str(first_file).endswith(".gz") else open
    try:
        with opener(first_file, "rt") as fh:
            header = fh.readline().strip().lstrip("@")
            parts = {k: v for p in header.split() if "=" in p for k, v in [p.split("=", 1)]}
            run_id      = parts.get("runid", hashlib.md5(str(fastq_dir.resolve()).encode()).hexdigest()[:12])
            flowcell_id = parts.get("flow_cell_id", "unknown")
            sample_id   = parts.get("sample_id", "unknown")
    except Exception:
        run_id = hashlib.md5(str(fastq_dir.resolve()).encode()).hexdigest()[:12]

    log(f"Run ID: {run_id} | Flowcell: {flowcell_id} | Sample: {sample_id}")

    conn = connect_db(db_path)
    if conn.execute("SELECT 1 FROM run WHERE run_id=?", (run_id,)).fetchone():
        log(f"Run {run_id} already ingested. Delete it first to re-ingest.")
        conn.close()
        return run_id

    now = datetime.utcnow().isoformat()

    # Group barcodes by target
    groups = group_barcodes_by_target(barcodes)
    log(f"Found {len(groups)} target(s): {', '.join(sorted(groups.keys()))}")

    with conn:
        conn.execute(
            """INSERT INTO run
               (run_id, flowcell_id, sample_id, folder_path,
                start_anchor, end_anchor, min_aa_len,
                mm_min_seq_id, mm_coverage, mm_cov_mode, use_linclust,
                n_targets, ingested_at)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (run_id, flowcell_id, sample_id, str(fastq_dir.resolve()),
             START, END, LENGTH,
             mm_min_seq_id, mm_coverage, mm_cov_mode, int(use_linclust),
             len(groups), now)
        )

    # Process each target
    for target, grp in sorted(groups.items()):
        log(f"\nProcessing target: {target}")
        tg1_info = grp["tg1"]
        rounds   = grp["rounds"]

        control_dir = Path(tg1_info["folder_path"]) if tg1_info else (ext_tg1_dir if ext_tg1_dir else None)
        onex_dir    = Path(rounds[1]["folder_path"]) if 1 in rounds else None
        twox_dir    = Path(rounds[2]["folder_path"]) if 2 in rounds else None

        if not onex_dir and not twox_dir:
            log(f"  [{target}] No panning rounds found — skipping")
            continue

        try:
            ingest_target(
                conn=conn, run_id=run_id, target=target,
                control_dir=control_dir, onex_dir=onex_dir, twox_dir=twox_dir,
                START=START, END=END, LENGTH=LENGTH,
                mm_min_seq_id=mm_min_seq_id, mm_coverage=mm_coverage,
                mm_cov_mode=mm_cov_mode, use_linclust=use_linclust,
                drop_unclustered=drop_unclustered,
                progress_cb=progress_cb,
            )
        except Exception as e:
            log(f"  [{target}] ERROR: {e}")

    log(f"\n✅ Ingest complete. Run ID: {run_id}")
    conn.close()
    return run_id


# ═══════════════════════════════════════════════════════════════════════════════
# LOAD FROM DB
# ═══════════════════════════════════════════════════════════════════════════════

def load_count_matrix(conn: sqlite3.Connection, run_id: str, target: str) -> pd.DataFrame:
    counts = sql_df(conn,
        "SELECT cluster_head, library_id, raw_count, norm_count FROM cluster_counts WHERE run_id=? AND target=?",
        (run_id, target))
    if counts.empty:
        return pd.DataFrame()
    raw_wide = counts.pivot_table(
        index="cluster_head", columns="library_id", values="raw_count",
        fill_value=0, aggfunc="sum").reset_index()
    norm_wide = counts.pivot_table(
        index="cluster_head", columns="library_id", values="norm_count",
        fill_value=0, aggfunc="sum").reset_index()
    norm_wide.columns = [f"{c}_norm" if c != "cluster_head" else c for c in norm_wide.columns]
    result = raw_wide.merge(norm_wide, on="cluster_head", how="left")
    if "2xpanned_norm" in result.columns:
        result = result.sort_values("2xpanned_norm", ascending=False).reset_index(drop=True)
    elif "1xpanned_norm" in result.columns:
        result = result.sort_values("1xpanned_norm", ascending=False).reset_index(drop=True)
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# STREAMLIT PAGES
# ═══════════════════════════════════════════════════════════════════════════════

def page_overview(conn: sqlite3.Connection):
    st.header("Overview")
    runs = sql_df(conn,
        "SELECT run_id, flowcell_id, sample_id, n_targets, start_anchor, end_anchor, ingested_at FROM run ORDER BY ingested_at DESC")
    if runs.empty:
        st.info("No runs ingested yet. Go to the Ingest page.")
        return
    st.dataframe(runs, use_container_width=True)

    run_id = st.selectbox("Select run for target details",
        runs["run_id"].tolist(),
        format_func=lambda r: runs[runs["run_id"]==r]["sample_id"].iloc[0])
    targets = sql_df(conn,
        """SELECT target, total_reads_control, total_reads_1xpanned, total_reads_2xpanned,
                  kept_filtered, n_clusters FROM target_run WHERE run_id=? ORDER BY target""",
        (run_id,))
    if not targets.empty:
        st.dataframe(targets, use_container_width=True)


def page_enrichment(conn: sqlite3.Connection):
    st.header("Enrichment & Candidates")
    runs = sql_df(conn, "SELECT run_id, sample_id FROM run ORDER BY ingested_at DESC")
    if runs.empty:
        st.info("No runs ingested yet.")
        return

    run_id = st.selectbox("Run",
        runs["run_id"].tolist(),
        format_func=lambda r: runs[runs["run_id"]==r]["sample_id"].iloc[0])

    targets_df = sql_df(conn, "SELECT target FROM target_run WHERE run_id=? ORDER BY target", (run_id,))
    if targets_df.empty:
        st.warning("No targets found for this run.")
        return

    target = st.selectbox("Target", targets_df["target"].tolist())

    col1, col2 = st.columns(2)
    with col1:
        baseline_threshold = st.number_input("Baseline threshold", min_value=1,
                                              value=DEFAULT_BASELINE_THRESHOLD, step=1)
    with col2:
        enrichment_threshold = st.number_input("Enrichment threshold (fold)", min_value=0.0001,
                                                value=DEFAULT_ENRICHMENT_THRESHOLD, step=1.0)

    show_sequences = st.checkbox("Show top cluster amino-acid sequences", value=True)
    top_n = st.slider("How many top clusters to show", 1, 100, 30, disabled=not show_sequences)

    cache_key = f"cm_{run_id}_{target}"
    if cache_key not in st.session_state:
        with st.spinner("Loading from database..."):
            st.session_state[cache_key] = load_count_matrix(conn, run_id, target)
    count_matrix = st.session_state[cache_key]

    if count_matrix.empty:
        st.warning("No cluster counts found for this target.")
        return

    # Use normalized counts (Test12.py)
    count_matrix_for_plots = count_matrix.copy()
    for lib in ["control", "1xpanned", "2xpanned"]:
        norm_col = f"{lib}_norm"
        if norm_col in count_matrix_for_plots.columns:
            count_matrix_for_plots[lib] = pd.to_numeric(count_matrix_for_plots[norm_col], errors="coerce")

    # Count matrix preview — add fold_enrichment column
    cm_display = count_matrix.copy()
    if "2xpanned_norm" in cm_display.columns and "control_norm" in cm_display.columns:
        cm_display["fold_enrichment"] = cm_display.apply(
            lambda r: round(float(r["2xpanned_norm"]) / float(r["control_norm"]), 2)
            if (r["control_norm"] is not None and float(r["control_norm"]) > 0) else None,
            axis=1
        )
        cm_display = cm_display.sort_values("fold_enrichment", ascending=False, na_position="last")
    elif "1xpanned_norm" in cm_display.columns and "control_norm" in cm_display.columns:
        cm_display["fold_enrichment"] = cm_display.apply(
            lambda r: round(float(r["1xpanned_norm"]) / float(r["control_norm"]), 2)
            if (r["control_norm"] is not None and float(r["control_norm"]) > 0) else None,
            axis=1
        )
        cm_display = cm_display.sort_values("fold_enrichment", ascending=False, na_position="last")

    st.subheader(f"Cluster count matrix — {target}")
    preferred = ["cluster_head","control","1xpanned","2xpanned","control_norm","1xpanned_norm","2xpanned_norm","fold_enrichment"]
    display_cols = [c for c in preferred if c in cm_display.columns] +                    [c for c in cm_display.columns if c not in preferred]
    st.dataframe(cm_display[display_cols].head(30), use_container_width=True)
    st.download_button("Download count matrix CSV",
                       cm_display[display_cols].to_csv(index=False).encode(),
                       f"{run_id[:8]}_{target}_count_matrix.csv", "text/csv")


    # Top sequences (Test12.py)
    if show_sequences:
        st.subheader(f"Top {top_n} cluster amino-acid sequences")
        top_ids = count_matrix["cluster_head"].astype(str).head(top_n).tolist()
        placeholders = ",".join(["?"]*len(top_ids))
        seqs_df = sql_df(conn,
            f"SELECT cluster_head, aa_sequence, aa_length FROM cluster_sequences WHERE run_id=? AND target=? AND cluster_head IN ({placeholders})",
            (run_id, target, *top_ids))
        if not seqs_df.empty:
            seqs_df["_order"] = seqs_df["cluster_head"].map({ch: i for i, ch in enumerate(top_ids)})
            seqs_df = seqs_df.sort_values("_order").drop("_order", axis=1)
            st.dataframe(seqs_df, use_container_width=True)
            fasta_text = seq_df_to_fasta_text(seqs_df)
            if fasta_text.strip():
                clipboard_copy_button(fasta_text, label="Copy nanobody sequences (FASTA)",
                                      key=f"{run_id}_{target}_top{top_n}")
                with st.expander("FASTA text (manual copy fallback)"):
                    st.text_area("FASTA", value=fasta_text, height=220)

    st.divider()

    # Abundance distribution (Test12.py)
    st.subheader("Abundance distribution by library")
    cluster_counts_long = sql_df(conn,
        "SELECT library_id, cluster_head, raw_count as n FROM cluster_counts WHERE run_id=? AND target=?",
        (run_id, target))
    if not cluster_counts_long.empty:
        with tempfile.TemporaryDirectory() as tmp:
            fig_abund, _ = plot_library_abundance_distribution(
                cluster_counts_long,
                out_pdf=Path(tmp)/"a.pdf", out_tsv=Path(tmp)/"a.tsv")
        st.pyplot(fig_abund, clear_figure=False)

    st.divider()

    # AA seq map for hover
    seqs_all = sql_df(conn,
        "SELECT cluster_head, aa_sequence FROM cluster_sequences WHERE run_id=? AND target=?",
        (run_id, target))
    aa_seq_map = dict(zip(seqs_all["cluster_head"], seqs_all["aa_sequence"])) if not seqs_all.empty else {}

    # Determine condition libs (may only have 1xpanned if no R2)
    condition_libs = [c for c in ["1xpanned", "2xpanned"] if c in count_matrix_for_plots.columns and count_matrix_for_plots[c].sum() > 0]

    # Interactive Plotly scatter (Test12.py)
    st.subheader("Differential abundance")
    diff_fig_plotly = plot_abundance_vs_differential_plotly(
        count_matrix=count_matrix_for_plots,
        baseline_lib="control",
        condition_libs=condition_libs,
        baseline_threshold=int(baseline_threshold),
        enrichment_threshold=float(enrichment_threshold),
        aa_seq_map=aa_seq_map,
    )
    st.plotly_chart(diff_fig_plotly, use_container_width=True)

    # Enriched candidates table
    with st.expander("Enriched candidates table + download"):
        with tempfile.TemporaryDirectory() as tmp:
            _, plot_data_diff = plot_abundance_vs_differential(
                count_matrix=count_matrix_for_plots,
                baseline_lib="control",
                condition_libs=condition_libs,
                baseline_threshold=int(baseline_threshold),
                enrichment_threshold=float(enrichment_threshold),
                out_pdf=Path(tmp)/"d.pdf", out_tsv=Path(tmp)/"d.tsv",
            )
        if not plot_data_diff.empty:
            enriched = plot_data_diff[plot_data_diff["highlight"]].copy()
            if not enriched.empty:
                st.write(f"**{len(enriched)} enriched clusters** (baseline ≥ {baseline_threshold}, fold ≥ {enrichment_threshold}×)")
                show_cols = [c for c in ["cluster_head","condition2","count1","count2","ratio"] if c in enriched.columns]
                enriched = enriched[show_cols].sort_values("ratio", ascending=False)
                if not seqs_all.empty:
                    enriched = enriched.merge(seqs_all, on="cluster_head", how="left")
                st.dataframe(enriched, use_container_width=True)
                st.download_button("Download enriched candidates CSV",
                                   enriched.to_csv(index=False).encode(),
                                   f"{run_id[:8]}_{target}_enriched.csv", "text/csv")
            else:
                st.info(f"No clusters meet thresholds (baseline ≥ {baseline_threshold}, fold ≥ {enrichment_threshold}×).")


def page_ingest(conn: sqlite3.Connection, db_path: Path):
    st.header("Ingest Run")
    st.caption("Point to a run folder — the app auto-detects all targets and processes each one.")

    folder_path = st.text_input("Path to run folder",
                                 placeholder="/Users/ishanghosh/Downloads/run5")

    if folder_path.strip():
        p = Path(folder_path.strip()).expanduser()
        if p.exists():
            barcodes = discover_barcodes(p)
            if barcodes:
                groups = group_barcodes_by_target(barcodes)
                st.success(f"Found {len(barcodes)} barcode folders → {len(groups)} target(s)")
                preview_rows = []
                for target, grp in sorted(groups.items()):
                    tg1 = grp["tg1"]["label"] if grp["tg1"] else "— (no TG1)"
                    r1  = grp["rounds"].get(1, {}).get("label", "—")
                    r2  = grp["rounds"].get(2, {}).get("label", "—")
                    preview_rows.append({"target": target, "TG1 (control)": tg1, "R1 (1xpanned)": r1, "R2 (2xpanned)": r2})
                st.dataframe(pd.DataFrame(preview_rows), use_container_width=True)
            else:
                st.error("No barcode folders found.")
        else:
            st.error("Path does not exist.")

    st.divider()
    START  = st.text_input("START anchor", value=DEFAULT_START)
    END    = st.text_input("END anchor",   value=DEFAULT_END)
    LENGTH = st.number_input("Minimum amino-acid length", min_value=1, value=DEFAULT_LENGTH, step=1)
    drop_unclustered = st.checkbox("Drop reads not in cluster TSV (recommended)", value=True)

    with st.expander("MMseqs2 parameters"):
        use_linclust  = st.checkbox("Use easy-linclust (fast) instead of easy-cluster",
                                    value=False,
                                    help="easy-cluster matches Test12.py exactly. easy-linclust is ~10x faster.")
        mm_min_seq_id = st.number_input("min_seq_id", 0.0, 1.0, 0.90, 0.01, format="%.2f")
        mm_coverage   = st.number_input("coverage",   0.0, 1.0, 0.90, 0.01, format="%.2f")
        mm_cov_mode   = st.number_input("cov_mode",   0, 5, 0, 1)

    # External TG1 selector
    st.divider()
    st.subheader("External TG1 (optional)")
    st.caption("Use this if the run has no TG1 barcode and TG1 comes from a different run.")

    ext_tg1_dir = None
    ext_tg1_source = st.radio("TG1 source", ["None", "Already ingested run", "Folder path"],
                               horizontal=True, key="ext_tg1_source")

    if ext_tg1_source == "Already ingested run":
        other_runs = sql_df(conn, "SELECT run_id, sample_id FROM run ORDER BY ingested_at DESC")
        if not other_runs.empty:
            sel_run = st.selectbox("Source run",
                other_runs["run_id"].tolist(),
                format_func=lambda r: other_runs[other_runs["run_id"]==r]["sample_id"].iloc[0],
                key="ext_tg1_run")
            tg1_paths = sql_df(conn,
                "SELECT DISTINCT control_path FROM target_run WHERE run_id=? AND control_path IS NOT NULL AND LOWER(control_path) LIKE '%tg1%'",
                (sel_run,))
            if not tg1_paths.empty:
                path_options = tg1_paths["control_path"].tolist()
                path_labels = [Path(p).name for p in path_options]
                sel_idx = st.selectbox("TG1 barcode folder",
                    range(len(path_labels)),
                    format_func=lambda i: path_labels[i],
                    key="ext_tg1_bc_sel")
                ext_tg1_dir = Path(path_options[sel_idx])
                st.caption(f"Path: `{ext_tg1_dir}`")
            else:
                st.info("No TG1 barcodes found in that run.")
        else:
            st.info("No runs ingested yet.")

    elif ext_tg1_source == "Folder path":
        tg1_path_input = st.text_input("TG1 barcode folder path",
            placeholder="/Users/ishanghosh/Downloads/run5/barcode06_P-aer_llama_TG1",
            key="ext_tg1_path")
        if tg1_path_input.strip():
            resolved = Path(tg1_path_input.strip()).expanduser()
            if resolved.exists():
                st.success(f"Found: `{resolved}`")
                ext_tg1_dir = resolved
            else:
                st.error("Path does not exist.")

    if ext_tg1_dir:
        st.info(f"External TG1 will be used for targets with no TG1 barcode: `{ext_tg1_dir.name}`")

    st.divider()
    if st.button("Start Ingest", type="primary", disabled=not folder_path.strip()):
        p = Path(folder_path.strip()).expanduser()
        if not p.exists():
            st.error("Path does not exist.")
            return
        log_area = st.empty()
        msgs = []
        def log(m):
            msgs.append(m)
            log_area.code("\n".join(msgs))
        try:
            run_id = ingest_run(
                fastq_dir=p, db_path=db_path,
                START=START, END=END, LENGTH=int(LENGTH),
                mm_min_seq_id=float(mm_min_seq_id),
                mm_coverage=float(mm_coverage),
                mm_cov_mode=int(mm_cov_mode),
                use_linclust=bool(use_linclust),
                drop_unclustered=bool(drop_unclustered),
                ext_tg1_dir=ext_tg1_dir,
                progress_cb=log,
            )
            st.success(f"✅ Ingest complete! Run ID: `{run_id}`")
        except Exception as e:
            st.error(f"Ingest failed: {e}")

    st.divider()
    st.subheader("Delete Run")
    runs = sql_df(conn, "SELECT run_id, sample_id FROM run ORDER BY ingested_at DESC")
    if not runs.empty:
        del_run = st.selectbox("Select run to delete",
            runs["run_id"].tolist(),
            format_func=lambda r: runs[runs["run_id"]==r]["sample_id"].iloc[0],
            key="del_run")
        if st.button("Delete run", type="secondary"):
            with conn:
                conn.execute("PRAGMA foreign_keys=OFF;")
                for tbl in ["cluster_sequences","cluster_counts","target_run","run"]:
                    conn.execute(f"DELETE FROM {tbl} WHERE run_id=?", (del_run,))
                conn.execute("PRAGMA foreign_keys=ON;")
            for k in list(st.session_state.keys()):
                if del_run in k:
                    del st.session_state[k]
            st.success("Deleted.")
            st.rerun()


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    st.set_page_config(page_title="MinION Nanobody Analysis v2", layout="wide")
    st.title("MinION Nanobody Analysis v2")

    db_path = Path(DEFAULT_DB).expanduser()
    conn = connect_db(db_path)

    page = st.sidebar.radio("Page", ["Overview", "Enrichment", "Ingest"])
    if page == "Overview":
        page_overview(conn)
    elif page == "Enrichment":
        page_enrichment(conn)
    elif page == "Ingest":
        page_ingest(conn, db_path)


if __name__ == "__main__":
    main()
