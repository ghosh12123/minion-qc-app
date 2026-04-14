"""
MinION Nanobody Analysis
"""

import gzip
import hashlib
import io
import multiprocessing as mp
import os
import json
import re
import shutil
import sqlite3
import subprocess
import tempfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components
from io import StringIO
import html as html_lib
import time
import requests
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
import matplotlib.colors as mcolors
from plotly.subplots import make_subplots

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

def build_cluster_counts_from_map(
    read_id_to_barcode: Dict[str, str],
    cluster_tsv: Path,
    drop_unclustered: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Optimized version of build_cluster_counts_and_matrix.
    Uses a pre-built {read_id: library_id} map instead of re-reading FASTQs.
    Produces identical output to the original function.
    """
    if not cluster_tsv.exists():
        raise FileNotFoundError(f"Cluster TSV not found: {cluster_tsv}")

    cluster_all = pd.read_csv(
        cluster_tsv, sep="\t", header=None,
        names=["cluster_head", "read_id"], dtype=str,
    )
    read_to_cluster = dict(
        zip(cluster_all["read_id"].astype(str), cluster_all["cluster_head"].astype(str))
    )

    counts = Counter()
    for read_id, lib_id in read_id_to_barcode.items():
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
        return cluster_counts, pd.DataFrame({"cluster_head": []})

    count_matrix = (
        cluster_counts.pivot_table(
            index="cluster_head", columns="library_id",
            values="n", fill_value=0, aggfunc="sum",
        ).reset_index()
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

        if float(baseline_threshold) > 0:
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


# Plotly version for interactive hover (ID + AA sequence)
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

    x_vline = np.log10(float(baseline_threshold)) if float(baseline_threshold) > 0 else None
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
        add_scatter(d_red, "darkred", 9, "sig in 2xpanned (red)", showlegend=(i == len(condition_libs)))

        if x_vline is not None:
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


# ClustalW (EMBL-EBI server) + MSA sorting + colored display
# -----------------------------
EBI_REST_ROOT = "https://www.ebi.ac.uk/Tools/services/rest"
CLUSTALW_SERVICE = "clustalw2"

AA20 = list("ACDEFGHIKLMNPQRSTVWY")
AA_CATEGORIES = AA20 + ["B", "Z", "J", "X", "U", "O", "*", "-"]


def _normalize_residue(ch: str) -> str:
    ch = (ch or "").upper()
    if ch == ".":
        ch = "-"
    if ch in AA_CATEGORIES:
        return ch
    if ch == "-":
        return "-"
    return "X"


def build_aa_color_map() -> dict[str, str]:
    """Deterministic per-amino-acid colors (20 distinct + a few extras)."""
    cmap = plt.get_cmap("tab20")
    colors = [mcolors.to_hex(cmap(i)) for i in range(20)]
    aa_colors = {aa: colors[i] for i, aa in enumerate(AA20)}
    aa_colors.update(
        {
            "-": "#FFFFFF",
            "X": "#BDBDBD",
            "B": "#CFCFCF",
            "Z": "#CFCFCF",
            "J": "#CFCFCF",
            "U": "#CFCFCF",
            "O": "#CFCFCF",
            "*": "#000000",
        }
    )
    return aa_colors


def _ebi_get(service: str, path: str, timeout_s: int = 60) -> str:
    url = f"{EBI_REST_ROOT}/{service}/{path.lstrip('/')}"
    r = requests.get(url, timeout=timeout_s)
    r.raise_for_status()
    return r.text


def _ebi_post_run(service: str, params: dict[str, str], timeout_s: int = 60) -> str:
    url = f"{EBI_REST_ROOT}/{service}/run"

    # IMPORTANT: do NOT send an "Expect" header; EBI may return 417 when it is present.
    headers = {
        "Accept": "text/plain",
        "User-Agent": "nanobody-seq-analyzer/1.0",
    }

    r = requests.post(url, data=params, headers=headers, timeout=timeout_s)

    # Helpful error message if EBI returns HTML/text describing the issue
    try:
        r.raise_for_status()
    except requests.HTTPError as e:
        raise requests.HTTPError(f"{e}\n\nEBI response body (first 1000 chars):\n{r.text[:1000]}") from e

    job_id = r.text.strip()
    if not job_id:
        raise RuntimeError(f"Empty job id returned by {url}\nResponse:\n{r.text[:1000]}")
    return job_id


def _ebi_poll_status(
    service: str, job_id: str, timeout_s: int = 900, poll_s: float = 2.0
) -> str:
    deadline = time.time() + float(timeout_s)
    last = ""
    while time.time() < deadline:
        last = _ebi_get(service, f"status/{job_id}", timeout_s=60).strip().upper()
        if last in {"FINISHED", "ERROR", "FAILURE", "NOT_FOUND"}:
            return last
        time.sleep(float(poll_s))
    raise TimeoutError(f"EBI job timed out: {service} job_id={job_id} last_status={last}")


def _ebi_result_types(service: str, job_id: str) -> list[str]:
    xml = _ebi_get(service, f"resulttypes/{job_id}", timeout_s=60)
    return re.findall(r"<identifier>(.*?)</identifier>", xml)


def _pick_result_type(result_types: list[str], *, prefer: list[str]) -> str | None:
    """
    prefer: list of substrings that must all be present (case-insensitive).
    """
    rt_lower = [(rt, rt.lower()) for rt in result_types]
    for rt, rtl in rt_lower:
        if all(p.lower() in rtl for p in prefer):
            return rt
    return None


@st.cache_data(show_spinner=False, ttl=7 * 24 * 3600)
def clustalw2_align_via_ebi(
    fasta_text: str,
    email: str,
    stype: str = "protein",
    timeout_s: int = 900,
) -> dict[str, str]:
    """
    Submit sequences to EMBL-EBI ClustalW2 server and return alignment text.
    Returns dict with keys:
      job_id, aln_fasta, aln_clustal (best-effort), result_types (json)
    """
    if not fasta_text or not fasta_text.strip():
        raise ValueError("Empty FASTA input to ClustalW2.")
    if not email or "@" not in email:
        raise ValueError("A valid email address is required by the EMBL-EBI ClustalW2 service.")

    job_id = _ebi_post_run(
        CLUSTALW_SERVICE,
        params={"email": email.strip(), "sequence": fasta_text, "stype": stype},
        timeout_s=60,
    )

    status = _ebi_poll_status(CLUSTALW_SERVICE, job_id, timeout_s=timeout_s, poll_s=2.0)
    if status != "FINISHED":
        raise RuntimeError(f"ClustalW2 job failed: status={status}, job_id={job_id}")

    rtypes = _ebi_result_types(CLUSTALW_SERVICE, job_id)

    fasta_type = next(
        (rt for rt in rtypes if "aln" in rt.lower() and "fasta" in rt.lower()),
        None,
    )

    clustal_type = next(
        (rt for rt in rtypes if "aln" in rt.lower() and "clustal" in rt.lower()),
        None,
    )
    if clustal_type is None:
        clustal_type = next((rt for rt in rtypes if rt.lower() == "aln"), None)

    aln_fasta = _ebi_get(CLUSTALW_SERVICE, f"result/{job_id}/{fasta_type}") if fasta_type else ""
    aln_clustal = _ebi_get(CLUSTALW_SERVICE, f"result/{job_id}/{clustal_type}") if clustal_type else ""

    if not (aln_fasta.strip() or aln_clustal.strip()):
        raise RuntimeError(f"EBI returned no alignment text. job_id={job_id} result_types={rtypes}")

    return {
        "job_id": job_id,
        "aln_fasta": aln_fasta,
        "aln_clustal": aln_clustal,
        "result_types": json.dumps(rtypes),
    }


def _parse_clustal_loose(text: str) -> MultipleSeqAlignment:
    """
    Tolerant CLUSTAL parser:
    - ignores consensus lines and non-sequence lines
    - accepts optional trailing position numbers
    """
    seq_chunks: dict[str, list[str]] = {}
    seen_header = False

    for raw in (text or "").splitlines():
        line = raw.rstrip("\n")

        if not seen_header:
            if line.lstrip().upper().startswith("CLUSTAL"):
                seen_header = True
            continue

        if not line.strip():
            continue

        # consensus lines start with whitespace in CLUSTAL
        if line[0].isspace():
            continue

        parts = line.split()
        if len(parts) < 2:
            continue

        sid = parts[0]
        chunk = parts[1]

        # skip junk lines
        if not chunk or set(chunk) <= set("*:."):
            continue

        seq_chunks.setdefault(sid, []).append(chunk)

    if not seq_chunks:
        raise ValueError("No sequence lines found in CLUSTAL text (loose parser).")

    records = [
        SeqRecord(Seq("".join(chunks)), id=sid, description="")
        for sid, chunks in seq_chunks.items()
    ]

    lengths = {len(r.seq) for r in records}
    if len(lengths) != 1:
        # Don’t silently “pad” this; it usually means the server output was truncated/odd.
        lens_preview = {r.id: len(r.seq) for r in records}
        raise ValueError(
            "CLUSTAL parse produced sequences of unequal lengths; output may be truncated or nonstandard.\n"
            f"Lengths: {lens_preview}"
        )

    return MultipleSeqAlignment(records)


def parse_msa_auto(aln_text: str) -> MultipleSeqAlignment:
    s = (aln_text or "").lstrip("\ufeff")  # strip BOM if present
    if not s.strip():
        raise ValueError("Empty alignment text returned by EBI.")

    head = s.lstrip()[:500].lower()
    if head.startswith("<!doctype html") or "<html" in head:
        raise ValueError("EBI returned HTML, not an alignment. First 400 chars:\n" + s[:400])

    looks_fasta = re.search(r"^>", s, flags=re.M) is not None
    looks_clustal = s.lstrip().upper().startswith("CLUSTAL")

    errors: dict[str, Exception] = {}

    # Prefer what it looks like, but try both
    fmt_order = []
    if looks_fasta:
        fmt_order.append("fasta")
    if looks_clustal:
        fmt_order.append("clustal")
    for fmt in ("fasta", "clustal"):
        if fmt not in fmt_order:
            fmt_order.append(fmt)

    # Try Biopython parsers first (use parse() not read() to tolerate extra sections)
    for fmt in fmt_order:
        try:
            alns = list(AlignIO.parse(StringIO(s), fmt))
            if alns:
                return alns[0]
        except Exception as e:
            errors[fmt] = e

    # If it’s CLUSTAL-like, fall back to tolerant parser
    if looks_clustal:
        try:
            return _parse_clustal_loose(s)
        except Exception as e:
            errors["clustal_loose"] = e

    msg = "Could not parse alignment as FASTA or CLUSTAL.\n"
    for k, e in errors.items():
        msg += f"- {k}: {type(e).__name__}: {e}\n"
    msg += "\nFirst 400 chars:\n" + s.strip()[:400]
    raise ValueError(msg)


def pairwise_identity(a: str, b: str) -> float:
    """Percent identity excluding columns where either sequence has a gap."""
    matches = 0
    compared = 0
    for ca, cb in zip(a, b):
        if ca == "-" or cb == "-":
            continue
        compared += 1
        if ca == cb:
            matches += 1
    return (matches / compared) if compared else 0.0

def sort_msa_by_pairwise(
    aln: MultipleSeqAlignment,
    mode: str = "Mean pairwise identity (desc)",
) -> tuple[MultipleSeqAlignment, pd.DataFrame, pd.DataFrame]:
    """
    Sort an MSA and return:
      - msa_sorted: MultipleSeqAlignment
      - msa_rank_df: per-sequence ranking/metrics table
      - pid_sorted: pairwise identity matrix (percent), sorted to match msa_sorted
    mode must be one of:
      - "Mean pairwise identity (desc)"
      - "Identity to top-1 (desc)"
      - "Server order"
    """
    if aln is None or len(aln) == 0:
        empty_rank = pd.DataFrame(
            columns=[
                "rank",
                "id",
                "orig_index",
                "mean_pairwise_identity_pct",
                "identity_to_top1_pct",
            ]
        )
        return aln, empty_rank, pd.DataFrame()

    ids = [rec.id for rec in aln]
    seqs = [str(rec.seq).upper().replace(".", "-") for rec in aln]

    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(f"MSA sequences are not all the same length: {sorted(lengths)}")

    n = len(seqs)
    pid = np.zeros((n, n), dtype=float)

    for i in range(n):
        pid[i, i] = 100.0 * pairwise_identity(seqs[i], seqs[i])
        for j in range(i + 1, n):
            v = 100.0 * pairwise_identity(seqs[i], seqs[j])
            pid[i, j] = v
            pid[j, i] = v

    if n > 1:
        mean_pid = (pid.sum(axis=1) - np.diag(pid)) / float(n - 1)
    else:
        mean_pid = np.array([np.nan], dtype=float)

    pid_to_top1 = pid[:, 0].copy()  # identity to first sequence in current (server) order

    if mode == "Server order":
        order = list(range(n))
    elif mode == "Mean pairwise identity (desc)":
        order = sorted(range(n), key=lambda i: (-mean_pid[i], i))
    elif mode == "Identity to top-1 (desc)":
        order = sorted(range(n), key=lambda i: (-pid_to_top1[i], i))
    else:
        raise ValueError(f"Unknown MSA_SORT_MODE: {mode}")

    msa_sorted = MultipleSeqAlignment([aln[i] for i in order])

    msa_rank_df = pd.DataFrame(
        {
            "rank": np.arange(1, n + 1, dtype=int),
            "id": [ids[i] for i in order],
            "orig_index": [i + 1 for i in order],  # 1-based original position
            "mean_pairwise_identity_pct": [float(mean_pid[i]) for i in order],
            "identity_to_top1_pct": [float(pid_to_top1[i]) for i in order],
        }
    )

    pid_df = pd.DataFrame(pid, index=ids, columns=ids)
    pid_sorted = pid_df.iloc[order, order].copy()
    pid_sorted.index = [ids[i] for i in order]
    pid_sorted.columns = [ids[i] for i in order]

    return msa_sorted, msa_rank_df, pid_sorted


def msa_to_text(aln: MultipleSeqAlignment, fmt: str = "fasta") -> str:
    buf = StringIO()
    AlignIO.write(aln, buf, fmt)
    return buf.getvalue()


def msa_to_colored_html(
    aln: MultipleSeqAlignment,
    block_size: int | None = 80,   # set to None (or 0) for single-row-per-seq
    max_id_chars: int = 36,
    aa_colors: dict[str, str] | None = None,
) -> str:
    if aa_colors is None:
        aa_colors = build_aa_color_map()

    ids = [r.id[:max_id_chars] for r in aln]
    seqs = [str(r.seq).upper().replace(".", "-") for r in aln]
    if not seqs:
        return "<div>No alignment available.</div>"

    L = len(seqs[0])
    name_w = min(max((len(i) for i in ids), default=0), max_id_chars)

    def colored_segment(seg: str) -> str:
        return "".join(
            f'<span style="color:{aa_colors.get(_normalize_residue(ch), "#000")};">'
            f"{html_lib.escape(ch)}</span>"
            for ch in seg
        )

    lines: list[str] = []

    # --- single-row mode (no wrapping; horizontal scroll) ---
    if block_size is None or int(block_size) <= 0:
        for sid, s in zip(ids, seqs):
            sid_pad = sid.ljust(name_w)
            lines.append(f"{html_lib.escape(sid_pad)}  {colored_segment(s)}")

        pre = "\n".join(lines)
        return f"""
        <!doctype html>
        <html>
        <head><meta charset="utf-8"></head>
        <body>
          <div style="
            font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace;
            font-size: 12px;
            line-height: 1.25;
            white-space: pre;          /* do not wrap */
            overflow-x: auto;          /* horizontal scroll */
            overflow-y: auto;          /* vertical scroll */
            overflow-wrap: normal;
            word-break: normal;
            height: 620px;
            border: 1px solid #ddd;
            padding: 10px;
            ">
{pre}
          </div>
        </body>
        </html>
        """

    # --- Original block mode (wrapped into blocks) ---
    for start in range(0, L, int(block_size)):
        end = min(L, start + int(block_size))
        for sid, s in zip(ids, seqs):
            sid_pad = sid.ljust(name_w)
            seg = s[start:end]
            lines.append(f"{html_lib.escape(sid_pad)}  {colored_segment(seg)}")
        lines.append("")

    pre = "\n".join(lines)
    return f"""
    <!doctype html>
    <html>
    <head><meta charset="utf-8"></head>
    <body>
      <div style="
        font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace;
        font-size: 12px;
        line-height: 1.25;
        white-space: pre;
        overflow: auto;
        height: 620px;
        border: 1px solid #ddd;
        padding: 10px;
        ">
{pre}
      </div>
    </body>
    </html>
    """


def plot_msa_heatmap_plotly(
    aln: MultipleSeqAlignment,
    title: str = "Top-50 nanobody MSA (ClustalW server; sorted by pairwise identity)",
    aa_colors: dict[str, str] | None = None,
    sticky_ids: set | None = None,
) -> go.Figure:
    if aa_colors is None:
        aa_colors = build_aa_color_map()
    if sticky_ids is None:
        sticky_ids = set()

    ids = [r.id for r in aln]
    seqs = [str(r.seq).upper().replace(".", "-") for r in aln]
    if not seqs:
        return go.Figure()

    L = len(seqs[0])
    cats = AA_CATEGORIES
    cat_to_code = {c: i for i, c in enumerate(cats)}
    K = len(cats)

    Z = np.zeros((len(seqs), L), dtype=int)
    T = np.empty((len(seqs), L), dtype=object)

    # Mark sticky rows in Y axis labels
    display_ids = [f"⚠️ {id_}" if id_ in sticky_ids else id_ for id_ in ids]

    for i, s in enumerate(seqs):
        for j, ch in enumerate(s):
            chn = _normalize_residue(ch)
            Z[i, j] = cat_to_code.get(chn, cat_to_code["X"])
            T[i, j] = chn

    # Discrete colorscale
    colorscale: list[tuple[float, str]] = []
    for i, c in enumerate(cats):
        col = aa_colors.get(c, "#000000")
        colorscale.append((i / K, col))
        colorscale.append(((i + 1) / K, col))

    max_id_len = max((len(x) for x in display_ids), default=10)
    left_margin = min(400, 9 * max_id_len + 40)

    fig = go.Figure()

    # Main heatmap
    fig.add_trace(go.Heatmap(
        z=Z,
        x=list(range(1, L + 1)),
        y=display_ids,
        text=T,
        hovertemplate="ID: %{y}<br>Align pos: %{x}<br>AA: %{text}<extra></extra>",
        colorscale=colorscale,
        zmin=-0.5,
        zmax=K - 0.5,
        showscale=False,
    ))

    # Add sticky indicator bar on the right side
    if sticky_ids:
        sticky_flags = [1 if id_ in sticky_ids else 0 for id_ in ids]
        sticky_colors = ["#fce8e8" if f else "rgba(0,0,0,0)" for f in sticky_flags]
        fig.add_trace(go.Bar(
            x=[L + 2] * len(display_ids),
            y=display_ids,
            orientation="h",
            marker=dict(color=sticky_colors, line=dict(width=0)),
            width=0.8,
            showlegend=True,
            name="⚠️ Sticky sequence",
            hovertemplate="⚠️ Sticky: %{y}<extra></extra>",
        ))

    fig.update_yaxes(autorange="reversed", title="Nanobody (cluster_head)")
    fig.update_xaxes(title="Alignment position")
    fig.update_layout(
        title=title,
        height=max(520, 18 * len(ids) + 180),
        margin=dict(l=left_margin, r=20, t=70, b=60),
        barmode="overlay",
    )
    return fig


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

    # Second pass: for targets with no TG1, try prefix match or global TG1
    tg1_barcodes = [bc for bc in barcodes if bc.get("round") == 0]
    # A barcode named just "TG1" (target == "TG1") is a global TG1 for all targets
    global_tg1 = next((bc for bc in tg1_barcodes if bc.get("target", "").upper() == "TG1"), None)
    for target, grp in groups.items():
        if grp["tg1"] is None:
            # Try prefix match first
            for tg1_bc in tg1_barcodes:
                tg1_target = tg1_bc.get("target", "")
                if tg1_target.upper() != "TG1" and target.startswith(tg1_target):
                    grp["tg1"] = tg1_bc
                    break
            # Fall back to global TG1
            if grp["tg1"] is None and global_tg1:
                grp["tg1"] = global_tg1

    # Remove targets with no rounds (pure TG1 barcodes with no panning)
    groups = {t: g for t, g in groups.items() if g["rounds"]}

    return groups


# ═══════════════════════════════════════════════════════════════════════════════
# DATABASE
# ═══════════════════════════════════════════════════════════════════════════════

DEFAULT_DB                   = "~/minion_nanobody.db"
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

-- Pre-clustering unique AA sequences per barcode (for sticky sequence detection)
CREATE TABLE IF NOT EXISTS raw_sequences (
  run_id        TEXT,
  target        TEXT,
  barcode       TEXT,
  aa_sequence   TEXT,
  aa_length     INTEGER,
  read_count    INTEGER,
  PRIMARY KEY (run_id, target, barcode, aa_sequence)
);

-- Sticky sequence results cache (persists across restarts)
CREATE TABLE IF NOT EXISTS sticky_cache (
  run_id            TEXT,
  target            TEXT,
  similarity_thresh REAL,
  cached_at         TEXT,
  result_json       TEXT,
  PRIMARY KEY (run_id, target, similarity_thresh)
);

-- Indexes for performance
CREATE INDEX IF NOT EXISTS idx_raw_sequences_aa ON raw_sequences(aa_sequence);
CREATE INDEX IF NOT EXISTS idx_raw_sequences_run ON raw_sequences(run_id);
CREATE INDEX IF NOT EXISTS idx_cluster_counts_run_target ON cluster_counts(run_id, target);
CREATE INDEX IF NOT EXISTS idx_cluster_sequences_run_target ON cluster_sequences(run_id, target);
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
# INGEST — whole run, one target at a time 
# ═══════════════════════════════════════════════════════════════════════════════

def _compute_target(args: dict) -> dict:
    """
    Pure compute worker — runs the full pipeline for one target and returns
    all results as a dict. No database writes. Safe to run in a subprocess.
    """
    import gzip, shutil, tempfile
    from pathlib import Path
    from collections import defaultdict
    from Bio import SeqIO
    from Bio.Seq import Seq
    import pandas as pd

    run_id      = args["run_id"]
    target      = args["target"]
    control_dir = Path(args["control_dir"]) if args.get("control_dir") else None
    onex_dir    = Path(args["onex_dir"])    if args.get("onex_dir")    else None
    twox_dir    = Path(args["twox_dir"])    if args.get("twox_dir")    else None
    START       = args["START"]
    END         = args["END"]
    LENGTH      = args["LENGTH"]
    mm_min_seq_id  = args["mm_min_seq_id"]
    mm_coverage    = args["mm_coverage"]
    mm_cov_mode    = args["mm_cov_mode"]
    use_linclust   = args["use_linclust"]
    drop_unclustered = args["drop_unclustered"]

    logs = []
    def log(msg):
        logs.append(f"  [{target}] {msg}")

    try:
        control_files = find_fastq_files(control_dir) if control_dir else []
        one_x_files   = find_fastq_files(onex_dir)    if onex_dir   else []
        two_x_files   = find_fastq_files(twox_dir)    if twox_dir   else []

        if not control_files and one_x_files and two_x_files:
            log("No TG1 — using R1 as control, R2 as 1xpanned")
            control_dir = onex_dir; control_files = one_x_files
            onex_dir = twox_dir;    one_x_files   = two_x_files
            twox_dir = None;        two_x_files   = []
        elif not control_files and one_x_files:
            log("No TG1 and only R1 — skipping")
            return {"target": target, "skipped": True, "logs": logs}

        all_inputs = control_files + one_x_files + two_x_files
        if not all_inputs:
            log("No FASTQ files — skipping")
            return {"target": target, "skipped": True, "logs": logs}

        log("Counting reads...")
        total_reads_control  = count_fastq_records(control_files)
        total_reads_1xpanned = count_fastq_records(one_x_files)
        total_reads_2xpanned = count_fastq_records(two_x_files)
        log(f"control={total_reads_control:,}  1x={total_reads_1xpanned:,}  2x={total_reads_2xpanned:,}")

        with tempfile.TemporaryDirectory(prefix=f"nb_{target}_") as tmp:
            tmp_path       = Path(tmp)
            trimmed_fasta  = tmp_path / "trimmed_proteins.fasta"
            filtered_fasta = tmp_path / "filtered.fasta"
            cluster_prefix = tmp_path / "clusters"
            cluster_tmp    = tmp_path / "mmseqs_tmp"
            cluster_tsv    = tmp_path / "clusters_cluster.tsv"

            log("Translating per-barcode sequences...")
            barcode_raw_seqs = {}
            all_lib_paths = {}
            if control_dir: all_lib_paths["control"]   = control_dir
            if onex_dir:    all_lib_paths["1xpanned"]  = onex_dir
            if twox_dir:    all_lib_paths["2xpanned"]  = twox_dir

            per_barcode_fastas = {}
            read_id_to_barcode: Dict[str, str] = {}  # for optimized count matrix
            for lib_id, lib_path in all_lib_paths.items():
                seq_counts = {}
                bc_fasta = tmp_path / f"trimmed_{lib_id}.fasta"
                with open(bc_fasta, "w") as bc_out:
                    for fq in find_fastq_files(lib_path):
                        opener = gzip.open(fq, "rt") if str(fq).endswith(".gz") else open(fq, "rt")
                        with opener as handle:
                            for record in SeqIO.parse(handle, "fastq"):
                                seq_str = str(record.seq).upper()
                                s_idx = seq_str.find(START)
                                e_idx = seq_str.find(END)
                                if s_idx == -1 or e_idx == -1 or e_idx <= s_idx:
                                    continue
                                sub = seq_str[s_idx: e_idx + len(END)]
                                codon_len = (len(sub) // 3) * 3
                                prot = str(Seq(sub[:codon_len]).translate(to_stop=False))
                                if len(prot) < LENGTH or "*" in prot[:LENGTH]:
                                    continue
                                seq_counts[prot] = seq_counts.get(prot, 0) + 1
                                bc_out.write(f">{record.id}\n{prot}\n")
                                read_id_to_barcode[record.id] = lib_id  # capture mapping
                barcode_raw_seqs[lib_id] = seq_counts
                per_barcode_fastas[lib_id] = bc_fasta
                log(f"  {lib_id}: {sum(seq_counts.values()):,} passing, {len(seq_counts):,} unique")

            log("Combining translated FASTAs...")
            with open(trimmed_fasta, "w") as combined_out:
                for bc_fasta in per_barcode_fastas.values():
                    with open(bc_fasta) as f_in:
                        shutil.copyfileobj(f_in, combined_out)

            log(f"Filtering (min AA length={LENGTH})...")
            kept6, discarded6 = filter_aa_fasta(trimmed_fasta, filtered_fasta, LENGTH=int(LENGTH))
            kept5 = kept6 + discarded6
            log(f"Kept {kept6:,}, discarded {discarded6:,}")

            if kept6 <= 0:
                log("No sequences passed filtering — skipping.")
                return {"target": target, "skipped": True, "logs": logs}

            raw_seq_rows = []
            for lib_id, seq_counts in barcode_raw_seqs.items():
                for aa_seq, count in seq_counts.items():
                    raw_seq_rows.append((run_id, target, lib_id, aa_seq, len(aa_seq), count))

            mode_str = "easy-linclust" if use_linclust else "easy-cluster"
            log(f"Clustering with MMseqs2 {mode_str}...")
            run_mmseqs_easy_cluster(
                filtered_fasta, cluster_prefix, cluster_tmp,
                min_seq_id=mm_min_seq_id, coverage=mm_coverage,
                cov_mode=mm_cov_mode, use_linclust=use_linclust,
            )
            if not cluster_tsv.exists():
                log("Cluster TSV not found — skipping")
                return {"target": target, "skipped": True, "logs": logs}
            log("Clustering complete")

            log("Building count matrix...")
            cluster_counts_df, count_matrix = build_cluster_counts_from_map(
                read_id_to_barcode=read_id_to_barcode,
                cluster_tsv=cluster_tsv,
                drop_unclustered=drop_unclustered,
            )
            for col in ["control", "1xpanned", "2xpanned"]:
                if col not in count_matrix.columns:
                    count_matrix[col] = 0

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

            log("Fetching AA sequences...")
            all_heads = count_matrix["cluster_head"].astype(str).tolist()
            seq_map = fetch_fasta_sequences_by_id(filtered_fasta, all_heads)

            # Prepare DB rows
            count_rows = []
            for _, row in count_matrix.iterrows():
                ch = str(row["cluster_head"])
                for lib, norm_col in [("control","control_norm"),("1xpanned","1xpanned_norm"),("2xpanned","2xpanned_norm")]:
                    raw = int(row.get(lib, 0) or 0)
                    nv  = row.get(norm_col, raw)
                    norm = float(nv) if nv is not pd.NA and nv is not None else float(raw)
                    count_rows.append((run_id, target, ch, lib, raw, norm))

            seq_rows = [(run_id, target, ch, seq, len(seq)) for ch, seq in seq_map.items()]

            return {
                "run_id": run_id,
                "target": target,
                "skipped": False,
                "logs": logs,
                "control_dir": str(control_dir) if control_dir else None,
                "onex_dir":    str(onex_dir)    if onex_dir    else None,
                "twox_dir":    str(twox_dir)    if twox_dir    else None,
                "total_reads_control":  total_reads_control,
                "total_reads_1xpanned": total_reads_1xpanned,
                "total_reads_2xpanned": total_reads_2xpanned,
                "kept5": kept5,
                "kept6": kept6,
                "n_clusters": n_clusters,
                "count_rows": count_rows,
                "seq_rows": seq_rows,
                "raw_seq_rows": raw_seq_rows,
            }

    except Exception as e:
        import traceback
        logs.append(f"  [{target}] ERROR: {e}\n{traceback.format_exc()}")
        return {"target": target, "skipped": True, "error": str(e), "logs": logs}


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
    """Sequential wrapper around _compute_target — used as fallback."""
    def log(msg):
        if progress_cb:
            progress_cb(msg)

    result = _compute_target({
        "run_id": run_id, "target": target,
        "control_dir": str(control_dir) if control_dir else None,
        "onex_dir":    str(onex_dir)    if onex_dir    else None,
        "twox_dir":    str(twox_dir)    if twox_dir    else None,
        "START": START, "END": END, "LENGTH": LENGTH,
        "mm_min_seq_id": mm_min_seq_id, "mm_coverage": mm_coverage,
        "mm_cov_mode": mm_cov_mode, "use_linclust": use_linclust,
        "drop_unclustered": drop_unclustered,
    })

    for msg in result.get("logs", []):
        log(msg)

    if result.get("skipped"):
        return

    _write_target_to_db(conn, result)
    log(f"  [{target}] Done ✅")


def _write_target_to_db(conn: sqlite3.Connection, result: dict):
    """Write computed target results to the database."""
    run_id = result["run_id"]
    target = result["target"]
    with conn:
        conn.execute(
            """INSERT OR REPLACE INTO target_run
               (run_id, target, control_path, onex_path, twox_path,
                total_reads_control, total_reads_1xpanned, total_reads_2xpanned,
                kept_trimmed, kept_filtered, n_clusters)
               VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
            (run_id, target,
             result["control_dir"], result["onex_dir"], result["twox_dir"],
             result["total_reads_control"], result["total_reads_1xpanned"], result["total_reads_2xpanned"],
             result["kept5"], result["kept6"], result["n_clusters"])
        )
        conn.executemany(
            "INSERT OR REPLACE INTO cluster_counts (run_id,target,cluster_head,library_id,raw_count,norm_count) VALUES (?,?,?,?,?,?)",
            result["count_rows"]
        )
        conn.executemany(
            "INSERT OR IGNORE INTO cluster_sequences (run_id,target,cluster_head,aa_sequence,aa_length) VALUES (?,?,?,?,?)",
            result["seq_rows"]
        )
        conn.executemany(
            "INSERT OR REPLACE INTO raw_sequences (run_id,target,barcode,aa_sequence,aa_length,read_count) VALUES (?,?,?,?,?,?)",
            result["raw_seq_rows"]
        )



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
    ext_overrides: Optional[Dict[str, Dict[str, Optional[Path]]]] = None,
    skip_targets: List[str] = None,
    progress_cb=None,
) -> str:
    def log(msg):
        if progress_cb:
            progress_cb(msg)

    skip_targets = skip_targets or []
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

    # Build args list for each target
    target_args = []
    for target, grp in sorted(groups.items()):
        if target in skip_targets:
            log(f"\nSkipping target: {target} (excluded by user)")
            continue
        tg1_info = grp["tg1"]
        rounds   = grp["rounds"]
        target_overrides = (ext_overrides or {}).get(target, {})
        control_dir = Path(tg1_info["folder_path"]) if tg1_info else (
            target_overrides.get("tg1") or (ext_tg1_dir if ext_tg1_dir else None))
        onex_dir = Path(rounds[1]["folder_path"]) if 1 in rounds else target_overrides.get("r1")
        twox_dir = Path(rounds[2]["folder_path"]) if 2 in rounds else target_overrides.get("r2")
        if not onex_dir and not twox_dir:
            log(f"  [{target}] No panning rounds found — skipping")
            continue
        target_args.append({
            "run_id": run_id, "target": target,
            "control_dir": str(control_dir) if control_dir else None,
            "onex_dir":    str(onex_dir)    if onex_dir    else None,
            "twox_dir":    str(twox_dir)    if twox_dir    else None,
            "START": START, "END": END, "LENGTH": LENGTH,
            "mm_min_seq_id": mm_min_seq_id, "mm_coverage": mm_coverage,
            "mm_cov_mode": mm_cov_mode, "use_linclust": use_linclust,
            "drop_unclustered": drop_unclustered,
        })

    if not target_args:
        log("No targets to process.")
        conn.close()
        return run_id

    # Run targets in parallel — use min(n_targets, cpu_count-1) workers
    n_workers = min(len(target_args), max(1, (os.cpu_count() or 2) - 1))
    log(f"\nProcessing {len(target_args)} target(s) in parallel ({n_workers} workers)...")

    if n_workers == 1 or len(target_args) == 1:
        # Single target — run directly (no multiprocessing overhead)
        results = [_compute_target(args) for args in target_args]
    else:
        # Multiple targets — run in parallel
        with mp.Pool(processes=n_workers) as pool:
            results = pool.map(_compute_target, target_args)

    # Write results to DB sequentially (safe — no concurrent writes)
    log("\nAll targets computed. Writing to database...")
    for result in results:
        for msg in result.get("logs", []):
            log(msg)
        if result.get("skipped"):
            if result.get("error"):
                log(f"  [{result['target']}] ERROR: {result['error']}")
            continue
        try:
            _write_target_to_db(conn, result)
            log(f"  [{result['target']}] Saved ✅")
        except Exception as e:
            log(f"  [{result['target']}] DB write ERROR: {e}")

    log(f"\n✅ Ingest complete. Run ID: {run_id}")
    conn.close()
    return run_id


# ═══════════════════════════════════════════════════════════════════════════════
# LOAD FROM DB
# ═══════════════════════════════════════════════════════════════════════════════

def find_sticky_sequences(
    conn: sqlite3.Connection,
    run_id: str,
    target: str,
    cluster_heads: List[str],
    similarity_threshold: float = 1.0,
) -> Dict[str, List[dict]]:
    """
    Check if cluster head sequences appear in raw_sequences from ANY other run.
    Uses a single batched SQL query for exact match (default).
    Returns {cluster_head: [{"run_id":..., "run_name":..., "target":..., "barcode":..., "read_count":...}]}
    """
    if not cluster_heads:
        return {}

    # Get AA sequences for all cluster heads in one query
    placeholders = ",".join(["?"] * len(cluster_heads))
    seqs_df = sql_df(conn,
        f"SELECT cluster_head, aa_sequence FROM cluster_sequences "
        f"WHERE run_id=? AND target=? AND cluster_head IN ({placeholders})",
        (run_id, target, *cluster_heads))

    if seqs_df.empty:
        return {}

    aa_to_cluster: Dict[str, str] = dict(zip(seqs_df["aa_sequence"], seqs_df["cluster_head"]))
    all_seqs = [s for s in aa_to_cluster.keys() if s]

    if not all_seqs:
        return {}

    results: Dict[str, List[dict]] = {}

    if similarity_threshold >= 1.0:
        # BATCHED exact match — single query for all sequences
        seq_placeholders = ",".join(["?"] * len(all_seqs))
        matches_df = sql_df(conn,
            f"SELECT aa_sequence, run_id, target, barcode, read_count "
            f"FROM raw_sequences "
            f"WHERE aa_sequence IN ({seq_placeholders}) AND run_id != ?",
            (*all_seqs, run_id))

        if not matches_df.empty:
            # Fetch run names in one query
            other_run_ids = matches_df["run_id"].unique().tolist()
            run_id_placeholders = ",".join(["?"] * len(other_run_ids))
            run_names_df = sql_df(conn,
                f"SELECT run_id, sample_id FROM run WHERE run_id IN ({run_id_placeholders})",
                tuple(other_run_ids))
            run_names = dict(zip(run_names_df["run_id"], run_names_df["sample_id"])) if not run_names_df.empty else {}

            for aa_seq, grp in matches_df.groupby("aa_sequence"):
                ch = aa_to_cluster.get(aa_seq)
                if not ch:
                    continue
                hit_list = []
                for _, m in grp.iterrows():
                    hit_list.append({
                        "run_id": m["run_id"],
                        "run_name": run_names.get(m["run_id"], m["run_id"][:8]),
                        "target": m["target"],
                        "barcode": m["barcode"],
                        "read_count": int(m["read_count"]),
                    })
                results[ch] = hit_list
    else:
        # Approximate match — per-sequence with length pre-filter
        def seq_identity(a: str, b: str) -> float:
            m = sum(ca == cb for ca, cb in zip(a, b))
            return m / max(len(a), len(b))

        # Fetch run names once
        run_names_df = sql_df(conn, "SELECT run_id, sample_id FROM run WHERE run_id != ?", (run_id,))
        run_names = dict(zip(run_names_df["run_id"], run_names_df["sample_id"])) if not run_names_df.empty else {}

        for aa_seq, ch in aa_to_cluster.items():
            min_len = int(len(aa_seq) * similarity_threshold)
            max_len = int(len(aa_seq) / similarity_threshold) + 1
            candidates = sql_df(conn,
                "SELECT run_id, target, barcode, aa_sequence, read_count "
                "FROM raw_sequences WHERE aa_length BETWEEN ? AND ? AND run_id != ?",
                (min_len, max_len, run_id))
            if candidates.empty:
                continue
            mask = candidates["aa_sequence"].apply(lambda s: seq_identity(aa_seq, s) >= similarity_threshold)
            matches = candidates[mask]
            if not matches.empty:
                hit_list = []
                for _, m in matches.iterrows():
                    hit_list.append({
                        "run_id": m["run_id"],
                        "run_name": run_names.get(m["run_id"], m["run_id"][:8]),
                        "target": m["target"],
                        "barcode": m["barcode"],
                        "read_count": int(m["read_count"]),
                    })
                results[ch] = hit_list

    return results




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


def generate_enrichment_pdf(
    conn: sqlite3.Connection,
    run_id: str,
    target: str,
    cm_display: pd.DataFrame,
    display_cols: list,
    sticky_map: dict,
    baseline_threshold: float,
    enrichment_threshold: float,
    diff_fig,
    abund_fig=None,
) -> bytes:
    """Generate a PDF report for a target enrichment analysis."""
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.lib import colors
    from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                     TableStyle, PageBreak, Image as RLImage)
    from reportlab.lib.enums import TA_CENTER, TA_LEFT
    import io as _io
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use("Agg")

    buf = _io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=letter,
                            leftMargin=0.75*inch, rightMargin=0.75*inch,
                            topMargin=0.75*inch, bottomMargin=0.75*inch)
    styles = getSampleStyleSheet()
    story = []

    title_style = ParagraphStyle("title", parent=styles["Title"], fontSize=16, spaceAfter=6)
    heading_style = ParagraphStyle("heading", parent=styles["Heading2"], fontSize=12, spaceAfter=4)
    normal_style = styles["Normal"]
    small_style = ParagraphStyle("small", parent=styles["Normal"], fontSize=8)

    # ── Title
    run_info = conn.execute("SELECT sample_id FROM run WHERE run_id=?", (run_id,)).fetchone()
    sample_id = run_info[0] if run_info else run_id[:8]
    story.append(Paragraph(f"Nanobody Enrichment Report", title_style))
    story.append(Paragraph(f"Run: {sample_id} &nbsp;|&nbsp; Target: {target}", styles["Heading3"]))
    story.append(Paragraph(f"Generated: {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}", small_style))
    story.append(Spacer(1, 12))

    # ── Ingest parameters
    story.append(Paragraph("Ingest Parameters", heading_style))
    run_params = conn.execute(
        "SELECT start_anchor, end_anchor, min_aa_len, mm_min_seq_id, mm_coverage, use_linclust FROM run WHERE run_id=?",
        (run_id,)).fetchone()
    if run_params:
        cluster_mode = "easy-linclust" if run_params[5] else "easy-cluster"
        param_data = [
            ["Min AA Length", "Min Seq Identity", "Coverage", "Clustering Mode"],
            [str(run_params[2]), f"{run_params[3]:.0%}", f"{run_params[4]:.0%}", cluster_mode],
        ]
        param_table = Table(param_data, colWidths=[1.5*inch]*4)
        param_table.setStyle(TableStyle([
            ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#f0f0f0")),
            ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
            ("FONTSIZE", (0,0), (-1,-1), 9),
            ("GRID", (0,0), (-1,-1), 0.5, colors.grey),
            ("ALIGN", (0,0), (-1,-1), "CENTER"),
            ("TOPPADDING", (0,0), (-1,-1), 4),
            ("BOTTOMPADDING", (0,0), (-1,-1), 4),
        ]))
        story.append(param_table)
    story.append(Spacer(1, 10))

    # ── Analysis thresholds
    story.append(Paragraph("Analysis Thresholds", heading_style))
    thresh_data = [
        ["Baseline Threshold (control_norm)", "Enrichment Threshold (fold)"],
        [str(baseline_threshold), str(enrichment_threshold)],
    ]
    thresh_table = Table(thresh_data, colWidths=[3*inch, 3*inch])
    thresh_table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#f0f0f0")),
        ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
        ("FONTSIZE", (0,0), (-1,-1), 9),
        ("GRID", (0,0), (-1,-1), 0.5, colors.grey),
        ("ALIGN", (0,0), (-1,-1), "CENTER"),
        ("TOPPADDING", (0,0), (-1,-1), 4),
        ("BOTTOMPADDING", (0,0), (-1,-1), 4),
    ]))
    story.append(thresh_table)
    story.append(Spacer(1, 10))

    # ── Target QC stats
    story.append(Paragraph("Target Statistics", heading_style))
    tr = conn.execute(
        "SELECT total_reads_control, total_reads_1xpanned, total_reads_2xpanned, kept_filtered, n_clusters FROM target_run WHERE run_id=? AND target=?",
        (run_id, target)).fetchone()
    if tr:
        qc_data = [
            ["Control Reads", "R1 Reads", "R2 Reads", "Sequences (post-filter)", "Clusters"],
            [f"{tr[0]:,}", f"{tr[1]:,}", f"{tr[2]:,}", f"{tr[3]:,}", f"{tr[4]:,}"],
        ]
        qc_table = Table(qc_data, colWidths=[1.4*inch]*5)
        qc_table.setStyle(TableStyle([
            ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#f0f0f0")),
            ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
            ("FONTSIZE", (0,0), (-1,-1), 9),
            ("GRID", (0,0), (-1,-1), 0.5, colors.grey),
            ("ALIGN", (0,0), (-1,-1), "CENTER"),
            ("TOPPADDING", (0,0), (-1,-1), 4),
            ("BOTTOMPADDING", (0,0), (-1,-1), 4),
        ]))
        story.append(qc_table)
    story.append(Spacer(1, 10))

    # ── Sticky summary
    n_sticky = len(sticky_map)
    sticky_color = colors.HexColor("#fce8e8") if n_sticky > 0 else colors.HexColor("#d4edda")
    sticky_text = f"Sticky sequences: {n_sticky} flagged (appear in other runs at {100}% identity)"
    story.append(Paragraph(sticky_text, ParagraphStyle("sticky", parent=normal_style,
                                                         backColor=sticky_color, borderPad=4)))
    story.append(Spacer(1, 10))

    # ── Differential abundance plot
    story.append(Paragraph("Differential Abundance Plot", heading_style))
    if diff_fig is not None:
        try:
            img_buf = _io.BytesIO()
            diff_fig.savefig(img_buf, format="png", dpi=150, bbox_inches="tight")
            img_buf.seek(0)
            img = RLImage(img_buf, width=6.5*inch, height=4*inch)
            story.append(img)
        except Exception:
            story.append(Paragraph("(Plot could not be rendered)", small_style))
    story.append(Spacer(1, 10))

    # ── Abundance distribution plot
    story.append(Paragraph("Abundance Distribution by Library", heading_style))
    if abund_fig is not None:
        try:
            img_buf2 = _io.BytesIO()
            abund_fig.savefig(img_buf2, format="png", dpi=150, bbox_inches="tight")
            img_buf2.seek(0)
            img2 = RLImage(img_buf2, width=4*inch, height=4*inch)
            story.append(img2)
        except Exception:
            story.append(Paragraph("(Plot could not be rendered)", small_style))
    story.append(Spacer(1, 10))

    # ── Count matrix table (top 30)
    story.append(PageBreak())
    story.append(Paragraph(f"Top Clusters by Fold Enrichment", heading_style))

    top_df = cm_display[display_cols].head(30).copy()
    # Flag sticky
    top_df["sticky"] = top_df["cluster_head"].apply(lambda x: "YES" if x in sticky_map else "")

    # Build table
    col_headers = ["#", "cluster_head"] + [c for c in display_cols if c != "cluster_head"] + ["sticky"]
    table_data = [col_headers]
    for i, (_, row) in enumerate(top_df.iterrows(), 1):
        r = [str(i), row["cluster_head"][:20] + "..."] +             [str(row.get(c, "")) for c in display_cols if c != "cluster_head"] +             [row.get("sticky", "")]
        table_data.append(r)

    col_widths = [0.3*inch, 2.2*inch] + [0.9*inch] * (len(display_cols) - 1) + [0.5*inch]
    cm_table = Table(table_data, colWidths=col_widths, repeatRows=1)
    table_style = [
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#f0f0f0")),
        ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
        ("FONTSIZE", (0,0), (-1,-1), 7),
        ("GRID", (0,0), (-1,-1), 0.3, colors.grey),
        ("ALIGN", (0,0), (-1,-1), "CENTER"),
        ("ALIGN", (1,0), (1,-1), "LEFT"),
        ("TOPPADDING", (0,0), (-1,-1), 2),
        ("BOTTOMPADDING", (0,0), (-1,-1), 2),
    ]
    # Highlight sticky rows
    for i, (_, row) in enumerate(top_df.iterrows(), 1):
        if row.get("sticky") == "YES":
            table_style.append(("BACKGROUND", (0,i), (-1,i), colors.HexColor("#fce8e8")))
    cm_table.setStyle(TableStyle(table_style))
    story.append(cm_table)
    story.append(Spacer(1, 10))

    # ── Top sequences
    story.append(PageBreak())
    story.append(Paragraph("Top Cluster Amino Acid Sequences", heading_style))
    top_ids = cm_display["cluster_head"].astype(str).head(30).tolist()
    placeholders = ",".join(["?"]*len(top_ids))
    seqs = conn.execute(
        f"SELECT cluster_head, aa_sequence, aa_length FROM cluster_sequences WHERE run_id=? AND target=? AND cluster_head IN ({placeholders})",
        (run_id, target, *top_ids)
    ).fetchall()
    seq_map_local = {s[0]: (s[1], s[2]) for s in seqs}

    seq_data = [["#", "cluster_head", "AA length", "Sticky", "AA sequence (first 60 AA)"]]
    for i, ch in enumerate(top_ids, 1):
        aa_seq, aa_len = seq_map_local.get(ch, ("", 0))
        is_sticky = "YES" if ch in sticky_map else ""
        seq_data.append([str(i), ch[:20]+"...", str(aa_len), is_sticky, (aa_seq or "")[:60]+"..."])

    seq_table = Table(seq_data, colWidths=[0.3*inch, 1.8*inch, 0.6*inch, 0.5*inch, 3.3*inch], repeatRows=1)
    seq_style = [
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#f0f0f0")),
        ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
        ("FONTSIZE", (0,0), (-1,-1), 6.5),
        ("GRID", (0,0), (-1,-1), 0.3, colors.grey),
        ("ALIGN", (0,0), (3,-1), "CENTER"),
        ("ALIGN", (4,0), (4,-1), "LEFT"),
        ("TOPPADDING", (0,0), (-1,-1), 2),
        ("BOTTOMPADDING", (0,0), (-1,-1), 2),
    ]
    for i, ch in enumerate(top_ids, 1):
        if ch in sticky_map:
            seq_style.append(("BACKGROUND", (0,i), (-1,i), colors.HexColor("#fce8e8")))
    seq_table.setStyle(TableStyle(seq_style))
    story.append(seq_table)

    doc.build(story)
    buf.seek(0)
    return buf.read()


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

    # Show ingest parameters for this run
    run_params = sql_df(conn,
        "SELECT start_anchor, end_anchor, min_aa_len, mm_min_seq_id, mm_coverage, mm_cov_mode, use_linclust FROM run WHERE run_id=?",
        (run_id,))
    if not run_params.empty:
        p = run_params.iloc[0]
        cluster_mode = "easy-linclust" if p["use_linclust"] else "easy-cluster"
        with st.expander("Ingest parameters used for this run", expanded=False):
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Min AA length", p["min_aa_len"])
            with col2:
                st.metric("Min seq identity", f"{p['mm_min_seq_id']:.0%}")
            with col3:
                st.metric("Coverage", f"{p['mm_coverage']:.0%}")

    col1, col2 = st.columns(2)
    with col1:
        baseline_threshold = st.number_input("Baseline threshold", min_value=0,
                                              value=DEFAULT_BASELINE_THRESHOLD, step=1)
    with col2:
        enrichment_threshold = st.number_input("Enrichment threshold (fold)", min_value=0.0001,
                                                value=DEFAULT_ENRICHMENT_THRESHOLD, step=1.0)

    # PDF export placeholder — above all sliders
    pdf_placeholder = st.empty()

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

    # Use normalized counts
    count_matrix_for_plots = count_matrix.copy()
    for lib in ["control", "1xpanned", "2xpanned"]:
        norm_col = f"{lib}_norm"
        if norm_col in count_matrix_for_plots.columns:
            count_matrix_for_plots[lib] = pd.to_numeric(count_matrix_for_plots[norm_col], errors="coerce")

    # Determine condition libs early (needed for PDF export and plots)
    condition_libs = [c for c in ["1xpanned", "2xpanned"] if c in count_matrix_for_plots.columns and count_matrix_for_plots[c].sum() > 0]

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

    # Compute sticky sequences early so we can highlight both tables
    sticky_threshold = st.slider(
        "Sticky sequence similarity threshold",
        min_value=0.80, max_value=1.0, value=1.0, step=0.01,
        format="%.2f",
        help="1.0 = exact match only. Lower values catch near-identical sequences across runs.",
        key=f"sticky_thresh_{run_id}_{target}"
    )

    top_ids = cm_display["cluster_head"].astype(str).head(top_n).tolist()
    other_runs = sql_df(conn, "SELECT run_id FROM run WHERE run_id != ?", (run_id,))
    if not other_runs.empty:
        sticky_cache_key = f"sticky_{run_id}_{target}_{sticky_threshold}"
        cached = conn.execute(
            "SELECT result_json FROM sticky_cache WHERE run_id=? AND target=? AND similarity_thresh=?",
            (run_id, target, float(sticky_threshold))
        ).fetchone()
        if cached:
            import json as _json
            sticky_map = _json.loads(cached[0])
        else:
            with st.spinner("Checking for sticky sequences across other runs..."):
                sticky_map = find_sticky_sequences(
                    conn, run_id, target, top_ids,
                    similarity_threshold=float(sticky_threshold)
                )
            import json as _json
            with conn:
                conn.execute(
                    """INSERT OR REPLACE INTO sticky_cache
                       (run_id, target, similarity_thresh, cached_at, result_json)
                       VALUES (?,?,?,?,?)""",
                    (run_id, target, float(sticky_threshold),
                     datetime.utcnow().isoformat(),
                     _json.dumps(sticky_map))
                )
    else:
        sticky_map = {}

    # Sticky cluster heads set
    sticky_cluster_heads = set(sticky_map.keys())

    def highlight_sticky_cm(row):
        if row.get("cluster_head") in sticky_cluster_heads:
            return ["background-color: #fce8e8; color: #666666"] * len(row)
        return [""] * len(row)

    st.subheader(f"Cluster count matrix — {target}")
    if sticky_cluster_heads:
        st.warning(f"⚠️ {len(sticky_cluster_heads)} sticky sequence(s) found — highlighted in red")
    else:
        st.success("✅ No sticky sequences found in other runs")
    preferred = ["cluster_head","control_norm","1xpanned_norm","2xpanned_norm","fold_enrichment"]
    display_cols = [c for c in preferred if c in cm_display.columns]
    cm_show = cm_display[display_cols].head(top_n)
    if sticky_cluster_heads:
        st.dataframe(cm_show.style.apply(highlight_sticky_cm, axis=1), use_container_width=True)
    else:
        st.dataframe(cm_show, use_container_width=True)
    st.download_button("Download count matrix CSV",
                       cm_display[display_cols].to_csv(index=False).encode(),
                       f"{run_id[:8]}_{target}_count_matrix.csv", "text/csv")


    # Top sequences with sticky detection
    if show_sequences:
        st.subheader(f"Top {top_n} cluster amino-acid sequences")

        top_ids = cm_display["cluster_head"].astype(str).head(top_n).tolist()
        placeholders = ",".join(["?"]*len(top_ids))
        seqs_df = sql_df(conn,
            f"SELECT cluster_head, aa_sequence, aa_length FROM cluster_sequences WHERE run_id=? AND target=? AND cluster_head IN ({placeholders})",
            (run_id, target, *top_ids))

        if not seqs_df.empty:
            seqs_df["_order"] = seqs_df["cluster_head"].map({ch: i for i, ch in enumerate(top_ids)})
            seqs_df = seqs_df.sort_values("_order").drop("_order", axis=1)

            # sticky_map already computed above

            # Build display dataframe with sticky info
            display_rows = []
            for _, row in seqs_df.iterrows():
                ch = row["cluster_head"]
                hits = sticky_map.get(ch, [])
                is_sticky = len(hits) > 0
                if hits:
                    # Summarize: "run_name / target (barcode) x3 reads"
                    hit_strs = []
                    for h in hits:
                        hit_strs.append(f"{h['run_name']} / {h['target']} ({h['barcode']}) — {h['read_count']} reads")
                    appears_in = " | ".join(hit_strs)
                    n_appearances = len(hits)
                else:
                    appears_in = ""
                    n_appearances = 0

                display_rows.append({
                    "cluster_head": ch,
                    "aa_sequence": row["aa_sequence"],
                    "aa_length": row["aa_length"],
                    "sticky": "⚠️ YES" if is_sticky else "",
                    "appears_in": appears_in,
                    "n_other_runs": n_appearances,
                })

            display_df = pd.DataFrame(display_rows)

            # Highlight sticky rows red
            def highlight_sticky(row):
                if row.get("sticky") == "⚠️ YES":
                    return ["background-color: #fce8e8; color: #666666"] * len(row)
                return [""] * len(row)

            n_sticky = display_df["sticky"].str.contains("YES", na=False).sum()
            if n_sticky > 0:
                st.warning(f"⚠️ {n_sticky} sticky sequence(s) found — highlighted in red")
            else:
                st.success("✅ No sticky sequences found in other runs")

            st.dataframe(
                display_df.style.apply(highlight_sticky, axis=1),
                use_container_width=True
            )

            fasta_text = seq_df_to_fasta_text(seqs_df)
            if fasta_text.strip():
                clipboard_copy_button(fasta_text, label="Copy nanobody sequences (FASTA)",
                                      key=f"{run_id}_{target}_top{top_n}")
                with st.expander("FASTA text (manual copy fallback)"):
                    st.text_area("FASTA", value=fasta_text, height=220)


    st.divider()
    # MSA panel
    with st.expander("Multiple sequence alignment (ClustalW server)", expanded=False):
        st.caption("Submits sequences to EMBL-EBI ClustalW2. Requires internet and a valid email.")
        RUN_MSA = st.checkbox("Run MSA", value=False, key=f"run_msa_{run_id}_{target}")

        MSA_SEQ_SOURCE = st.radio(
            "Sequences to align",
            ["Top N by count (default)", "Enriched candidates only (above threshold)"],
            horizontal=True,
            disabled=not RUN_MSA,
            key=f"msa_source_{run_id}_{target}"
        )
        MSA_N = st.slider("Max sequences to align", 2, 100, 50,
                          disabled=not RUN_MSA or MSA_SEQ_SOURCE != "Top N by count (default)",
                          key=f"msa_n_{run_id}_{target}")

        CLUSTALW_EMAIL = st.text_input("Email for EMBL-EBI submission (required)",
                                        value="", disabled=not RUN_MSA,
                                        key=f"msa_email_{run_id}_{target}")
        MSA_SORT_MODE = st.selectbox("Sort sequences by",
            ["Mean pairwise identity (desc)", "Identity to top-1 (desc)", "Server order"],
            index=0, disabled=not RUN_MSA, key=f"msa_sort_{run_id}_{target}")

        if RUN_MSA:
            if not CLUSTALW_EMAIL.strip() or "@" not in CLUSTALW_EMAIL:
                st.warning("Enter a valid email address to run ClustalW2.")
            else:
                # Determine which sequences to align
                if MSA_SEQ_SOURCE == "Enriched candidates only (above threshold)":
                    # Use clusters that meet the enrichment threshold
                    if "fold_enrichment" in cm_display.columns:
                        enriched_ids = cm_display[
                            cm_display["fold_enrichment"].notna() &
                            (cm_display["fold_enrichment"] >= float(enrichment_threshold)) &
                            (cm_display["control_norm"] >= float(baseline_threshold))
                        ]["cluster_head"].astype(str).head(100).tolist()
                    else:
                        enriched_ids = []
                    if len(enriched_ids) < 2:
                        st.warning(f"Only {len(enriched_ids)} enriched candidate(s) found above thresholds — need at least 2. Try lowering the thresholds.")
                        msa_ids = []
                    else:
                        msa_ids = enriched_ids
                        st.caption(f"Aligning {len(msa_ids)} enriched candidates (fold ≥ {enrichment_threshold}×, baseline ≥ {baseline_threshold})")
                else:
                    msa_ids = cm_display["cluster_head"].astype(str).head(MSA_N).tolist()
                    st.caption(f"Aligning top {len(msa_ids)} sequences by count")

                if msa_ids:
                    placeholders = ",".join(["?"]*len(msa_ids))
                    seqs50 = sql_df(conn,
                        f"SELECT cluster_head, aa_sequence FROM cluster_sequences WHERE run_id=? AND target=? AND cluster_head IN ({placeholders})",
                        (run_id, target, *msa_ids))
                    if seqs50.empty or len(seqs50) < 2:
                        st.info("Need at least 2 sequences for alignment.")
                    else:
                        seqs50["_order"] = seqs50["cluster_head"].map({ch: i for i, ch in enumerate(msa_ids)})
                        seqs50 = seqs50.sort_values("_order").drop("_order", axis=1)
                        fasta_in = seq_df_to_fasta_text(seqs50)
                    try:
                        with st.spinner("Submitting to ClustalW2 server..."):
                            aln_result = clustalw2_align_via_ebi(
                                fasta_text=fasta_in,
                                email=CLUSTALW_EMAIL.strip(),
                                stype="protein",
                            )
                        aln_text = (aln_result.get("aln_fasta", "") or "").strip() or (aln_result.get("aln_clustal", "") or "").strip()
                        msa = parse_msa_auto(aln_text)
                        msa_sorted, msa_rank_df, pid_sorted = sort_msa_by_pairwise(msa, mode=MSA_SORT_MODE)

                        tabs = st.tabs(["Alignment viewer", "Pairwise identity", "Raw ClustalW output"])
                        with tabs[0]:
                            # Pass sticky_ids from the already-computed sticky_map
                            msa_sticky_ids = set(sticky_map.keys()) if sticky_map else set()
                            # Build a set of prefixes that match MSA truncated IDs
                            msa_display_ids = {str(r.id) for r in msa_sorted}
                            msa_sticky_matched = set()
                            for did in msa_display_ids:
                                for sid in msa_sticky_ids:
                                    if sid.startswith(did) or did.startswith(sid):
                                        msa_sticky_matched.add(did)
                                        break
                            fig_msa = plot_msa_heatmap_plotly(msa_sorted, sticky_ids=msa_sticky_matched)
                            st.plotly_chart(fig_msa, use_container_width=True)
                            with st.expander("Colored alignment text (scroll)", expanded=False):
                                msa_html = msa_to_colored_html(msa_sorted, block_size=None)
                                components.html(msa_html, height=660, scrolling=True)
                        with tabs[1]:
                            # Add sticky column to msa_rank_df
                            rank_display = msa_rank_df.copy()
                            # Match on prefix since ClustalW truncates IDs
                            def _is_sticky(id_val, sticky_set):
                                if id_val in sticky_set:
                                    return "⚠️ YES"
                                # Try prefix match (ClustalW truncates long IDs)
                                for sid in sticky_set:
                                    if sid.startswith(id_val) or id_val.startswith(sid):
                                        return "⚠️ YES"
                                return ""
                            rank_display["sticky"] = rank_display["id"].apply(
                                lambda x: _is_sticky(x, msa_sticky_ids)
                            )
                            # Reorder columns to put sticky first
                            cols = ["rank", "id", "sticky"] + [c for c in rank_display.columns if c not in ["rank","id","sticky"]]
                            rank_display = rank_display[cols]

                            def highlight_sticky_row(row):
                                if row.get("sticky") == "⚠️ YES":
                                    return ["background-color: #fce8e8; color: #666666"] * len(row)
                                return [""] * len(row)

                            st.dataframe(
                                rank_display.style.apply(highlight_sticky_row, axis=1),
                                use_container_width=True
                            )
                            with st.expander("Full pairwise identity matrix"):
                                st.dataframe(pid_sorted, use_container_width=True)
                        with tabs[2]:
                            raw = aln_result.get("aln_clustal", "").strip()
                            if raw:
                                st.code(raw, language="text")
                            else:
                                st.info("No CLUSTAL-formatted output returned (FASTA alignment was used).")
                    except Exception as e:
                        st.error(f"MSA failed: {e}")


    # Abundance distribution 
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

    # Interactive Plotly scatter 
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

    # Populate PDF export button now that all data is available
    with pdf_placeholder.container():
        if st.button("Export enrichment report (PDF)", key=f"pdf_{run_id}_{target}"):
            with st.spinner("Generating PDF..."):
                try:
                    with tempfile.TemporaryDirectory() as tmp:
                        diff_fig_for_pdf, _ = plot_abundance_vs_differential(
                            count_matrix=count_matrix_for_plots,
                            baseline_lib="control",
                            condition_libs=condition_libs,
                            baseline_threshold=int(baseline_threshold),
                            enrichment_threshold=float(enrichment_threshold),
                            out_pdf=Path(tmp)/"d.pdf",
                            out_tsv=Path(tmp)/"d.tsv",
                        )
                        # Generate abundance figure
                        cluster_counts_long = sql_df(conn,
                            "SELECT library_id, cluster_head, raw_count as n FROM cluster_counts WHERE run_id=? AND target=?",
                            (run_id, target))
                        if not cluster_counts_long.empty:
                            abund_fig_for_pdf, _ = plot_library_abundance_distribution(
                                cluster_counts_long,
                                out_pdf=Path(tmp)/"a.pdf",
                                out_tsv=Path(tmp)/"a.tsv"
                            )
                        else:
                            abund_fig_for_pdf = None
                    pdf_bytes = generate_enrichment_pdf(
                        conn=conn, run_id=run_id, target=target,
                        cm_display=cm_display, display_cols=display_cols,
                        sticky_map=sticky_map,
                        baseline_threshold=baseline_threshold,
                        enrichment_threshold=enrichment_threshold,
                        diff_fig=diff_fig_for_pdf,
                        abund_fig=abund_fig_for_pdf,
                    )
                    st.download_button(
                        "Download PDF",
                        pdf_bytes,
                        f"{run_id[:8]}_{target}_enrichment_report.pdf",
                        "application/pdf",
                        key=f"pdf_dl_{run_id}_{target}"
                    )
                except Exception as e:
                    st.error(f"PDF generation failed: {e}")

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
    if "skip_targets" not in st.session_state:
        st.session_state["skip_targets"] = []
    st.caption("Point to a run folder — the app auto-detects all targets and processes each one.")

    folder_path = st.text_input("Path to run folder",
                                 placeholder="/Users/username/Downloads/run5")

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

                # Target skip checkboxes
                st.caption("Uncheck any targets you want to skip:")
                skip_targets = []
                cols = st.columns(min(4, len(groups)))
                for idx, t in enumerate(sorted(groups.keys())):
                    with cols[idx % len(cols)]:
                        if not st.checkbox(t, value=True, key=f"ingest_target_{t}"):
                            skip_targets.append(t)
                if skip_targets:
                    st.warning(f"Will skip: {', '.join(skip_targets)}")
                else:
                    st.session_state["skip_targets"] = []
                st.session_state["skip_targets"] = skip_targets
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
                                    value=False)
        mm_min_seq_id = st.number_input("min_seq_id", 0.0, 1.0, 0.90, 0.01, format="%.2f")
        mm_coverage   = st.number_input("coverage",   0.0, 1.0, 0.90, 0.01, format="%.2f")
        mm_cov_mode   = st.number_input("cov_mode",   0, 5, 0, 1)

    # ── External library overrides ──────────────────────────────────────────────
    st.divider()
    st.subheader("External library overrides (optional)")
    st.caption(
        "For any missing TG1, R1, or R2 in a target group, you can supply a barcode folder "
        "from an already-ingested run or provide a direct folder path."
    )

    # Build a helper widget that returns a Path or None for a given slot
    def library_selector(label: str, key_prefix: str) -> Optional[Path]:
        source = st.radio(
            f"{label} source",
            ["None", "Already ingested run", "Folder path"],
            horizontal=True,
            key=f"{key_prefix}_source",
        )
        if source == "Already ingested run":
            other_runs = sql_df(conn, "SELECT run_id, sample_id FROM run ORDER BY ingested_at DESC")
            if other_runs.empty:
                st.info("No runs ingested yet.")
                return None
            sel_run = st.selectbox(
                "Source run",
                other_runs["run_id"].tolist(),
                format_func=lambda r: other_runs[other_runs["run_id"] == r]["sample_id"].iloc[0],
                key=f"{key_prefix}_run",
            )
            # Get all barcode folders stored for that run
            all_paths = sql_df(conn,
                """SELECT DISTINCT control_path as path FROM target_run WHERE run_id=? AND control_path IS NOT NULL
                   UNION
                   SELECT DISTINCT onex_path   as path FROM target_run WHERE run_id=? AND onex_path IS NOT NULL
                   UNION
                   SELECT DISTINCT twox_path   as path FROM target_run WHERE run_id=? AND twox_path IS NOT NULL""",
                (sel_run, sel_run, sel_run))
            if all_paths.empty:
                st.info("No barcode folders found for that run.")
                return None
            path_options = all_paths["path"].tolist()
            path_labels  = [Path(p).name for p in path_options]
            sel_idx = st.selectbox(
                "Barcode folder",
                range(len(path_labels)),
                format_func=lambda i: path_labels[i],
                key=f"{key_prefix}_bc",
            )
            chosen = Path(path_options[sel_idx])
            st.caption(f"`{chosen}`")
            return chosen
        elif source == "Folder path":
            raw = st.text_input(
                f"{label} folder path",
                placeholder="/Users/username/Downloads/run5/barcode06_P-aer_llama_TG1",
                key=f"{key_prefix}_path",
            )
            if raw.strip():
                resolved = Path(raw.strip()).expanduser()
                if resolved.exists():
                    st.success(f"Found: `{resolved}`")
                    return resolved
                else:
                    st.error("Path does not exist.")
        return None

    # Show override widgets only for targets that are missing one or more slots
    ext_overrides: Dict[str, Dict[str, Optional[Path]]] = {}  # {target: {slot: path}}

    if folder_path.strip():
        p_check = Path(folder_path.strip()).expanduser()
        if p_check.exists():
            barcodes_check = discover_barcodes(p_check)
            if barcodes_check:
                groups_check = group_barcodes_by_target(barcodes_check)
                missing_targets = {
                    t: g for t, g in groups_check.items()
                    if g["tg1"] is None or 1 not in g["rounds"] or 2 not in g["rounds"]
                }
                if missing_targets:
                    for target, grp in sorted(missing_targets.items()):
                        with st.expander(f"⚠️  {target} — missing libraries"):
                            ext_overrides[target] = {}
                            if grp["tg1"] is None:
                                st.markdown("**TG1 (control) — not found in run folder**")
                                ext_overrides[target]["tg1"] = library_selector("TG1", f"{target}_tg1")
                            if 1 not in grp["rounds"]:
                                st.markdown("**R1 (1xpanned) — not found in run folder**")
                                ext_overrides[target]["r1"] = library_selector("R1", f"{target}_r1")
                            if 2 not in grp["rounds"]:
                                st.markdown("**R2 (2xpanned) — not found in run folder**")
                                ext_overrides[target]["r2"] = library_selector("R2", f"{target}_r2")
                else:
                    st.success("All targets have complete TG1 / R1 / R2 — no overrides needed.")

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
                ext_tg1_dir=None,
                ext_overrides=ext_overrides,
                skip_targets=list(st.session_state.get("skip_targets", [])),
                progress_cb=log,
            )
            st.success(f"✅ Ingest complete! Run ID: `{run_id}`")
            # Clear sticky DB cache — new run affects results for all existing targets
            with conn:
                conn.execute("DELETE FROM sticky_cache")
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
            try:
                del_conn = connect_db(db_path)
                del_conn.execute("PRAGMA foreign_keys=OFF;")
                for tbl in ["cluster_sequences","cluster_counts","target_run",
                            "raw_sequences","sticky_cache","run"]:
                    del_conn.execute(f"DELETE FROM {tbl} WHERE run_id=?", (del_run,))
                del_conn.execute("DELETE FROM sticky_cache")
                del_conn.commit()
                del_conn.execute("PRAGMA foreign_keys=ON;")
                del_conn.close()
                for k in list(st.session_state.keys()):
                    if del_run in k:
                        del st.session_state[k]
                st.success("Deleted.")
                st.rerun()
            except Exception as e:
                st.error(f"Delete failed: {e}. Try restarting the app first.")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    st.set_page_config(page_title="MinION Nanobody Analysis", layout="wide")
    st.title("MinION Nanobody Analysis")

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
