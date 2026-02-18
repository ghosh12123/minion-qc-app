# app.py
# MinION QC + SQLite + Streamlit UI (OPTIMIZED VERSION)
#
# Performance improvements:
# - Single-pass BAM processing (combined stats + MAPQ)
# - Unified streaming with configurable filters
# - Parallel BAM processing with multiprocessing
# - Batch inserts for sticky sequences
# - Extracted config constants
#
# Requirements (in your venv):
#   pip install streamlit pandas
# System tools:
#   samtools must be on PATH
#
# Run:
#   streamlit run app.py

import os
import re
import io
import sys
import time
import json
import math
import shutil
import hashlib
import sqlite3
import tempfile
import datetime as dt
import subprocess
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Iterable, Iterator
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import streamlit as st
import zipfile
import psutil


import anthropic


# Anthropic API key for chatbot (hardcoded for all users)
ANTHROPIC_API_KEY = "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"


# -----------------------------
# System RAM Detection
# -----------------------------
def get_available_ram_gb() -> float:
    """Get available system RAM in GB."""
    return psutil.virtual_memory().available / (1024**3)


def get_total_ram_gb() -> float:
    """Get total system RAM in GB."""
    return psutil.virtual_memory().total / (1024**3)


def recommend_max_sequences() -> int:
    """
    Recommend max sequences to cluster based on available RAM.
    MMseqs2 needs ~0.5-1GB per 10K sequences for protein clustering.
    """
    available_gb = get_available_ram_gb()
    
    if available_gb >= 28:
        return 100_000  # 32GB+ system
    elif available_gb >= 12:
        return 50_000   # 16GB system
    elif available_gb >= 6:
        return 20_000   # 8GB system
    elif available_gb >= 3:
        return 10_000   # 4GB system
    else:
        return 5_000    # <4GB system (minimal)


# -----------------------------
# Config Constants
# -----------------------------
APP_TITLE = "MinION Run QC Database"
DEFAULT_DB = "minion_master.db"

# Database settings
DB_TIMEOUT_SECONDS = 60
DB_BUSY_TIMEOUT_MS = 60000

# Processing settings
DELETE_CHUNK_SIZE = 50_000
STICKY_BATCH_SIZE = 10_000
MAX_WORKERS = 4  # for parallel BAM processing

# ONT BAM flag constants (SAM)
FLAG_PAIRED = 0x1
FLAG_PROPER_PAIR = 0x2
FLAG_UNMAPPED = 0x4
FLAG_MATE_UNMAPPED = 0x8
FLAG_REVERSE = 0x10
FLAG_MATE_REVERSE = 0x20
FLAG_READ1 = 0x40
FLAG_READ2 = 0x80
FLAG_SECONDARY = 0x100
FLAG_QCFAIL = 0x200
FLAG_DUP = 0x400
FLAG_SUPPLEMENTARY = 0x800


# -----------------------------
# Small utilities
# -----------------------------
def utc_iso() -> str:
    return dt.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"


def na(x) -> str:
    if x is None:
        return "NA"
    s = str(x).strip()
    return s if s else "NA"


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def human_bytes(n: int) -> str:
    if n == 0:
        return "0 B"
    units = ["B", "KB", "MB", "GB", "TB"]
    k = 1024
    i = int(math.floor(math.log(n, k)))
    i = min(i, len(units) - 1)
    return f"{n / (k ** i):.2f} {units[i]}"


def run_cmd(cmd: List[str], cwd: Optional[str] = None, timeout: Optional[int] = None) -> subprocess.CompletedProcess:
    """
    Runs a command; raises RuntimeError on non-zero exit.
    Returns CompletedProcess with stdout/stderr as text.
    """
    p = subprocess.run(
        cmd,
        cwd=cwd,
        timeout=timeout,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if p.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            f"  {' '.join(cmd)}\n\n"
            f"Exit code: {p.returncode}\n\n"
            f"STDOUT:\n{p.stdout}\n\n"
            f"STDERR:\n{p.stderr}\n"
        )
    return p


def check_samtools() -> None:
    try:
        run_cmd(["samtools", "--version"], timeout=10)
    except Exception as e:
        st.error("samtools is required but not found/working on PATH.")
        st.code(str(e))
        st.stop()


# -----------------------------
# SQLite schema + helpers
# -----------------------------
SCHEMA_SQL = """
PRAGMA journal_mode=WAL;

CREATE TABLE IF NOT EXISTS run (
  run_id TEXT PRIMARY KEY,
  flowcell_id TEXT,
  device_id TEXT,
  device_type TEXT,
  exp_start_time TEXT,
  sequencing_kit TEXT,
  flowcell_type TEXT,
  experiment_type TEXT,
  sample_id TEXT,
  created_at TEXT
);

CREATE TABLE IF NOT EXISTS run_qc (
  run_id TEXT PRIMARY KEY,
  n_reads_total INTEGER,
  mean_read_length REAL,
  n_primary INTEGER,
  n_secondary INTEGER,
  n_supplementary INTEGER,
  mapped_reads_total INTEGER,
  mapping_pct REAL,
  created_at TEXT,
  FOREIGN KEY(run_id) REFERENCES run(run_id)
);

CREATE TABLE IF NOT EXISTS barcode_qc (
  run_id TEXT,
  barcode TEXT,
  total_reads INTEGER,
  mapped_reads_total INTEGER,
  mapping_pct REAL,
  mean_read_length REAL,
  mapq_mean REAL,
  created_at TEXT,
  PRIMARY KEY (run_id, barcode),
  FOREIGN KEY(run_id) REFERENCES run(run_id)
);

CREATE TABLE IF NOT EXISTS alignment_groups (
  run_id TEXT,
  ref_name TEXT,
  read_count INTEGER,
  mapq_mean REAL,
  created_at TEXT,
  PRIMARY KEY (run_id, ref_name),
  FOREIGN KEY(run_id) REFERENCES run(run_id)
);

CREATE TABLE IF NOT EXISTS alignment_groups_by_barcode (
  run_id TEXT,
  barcode TEXT,
  ref_name TEXT,
  read_count INTEGER,
  mapq_mean REAL,
  created_at TEXT,
  PRIMARY KEY (run_id, barcode, ref_name),
  FOREIGN KEY(run_id) REFERENCES run(run_id)
);

CREATE TABLE IF NOT EXISTS sticky_sequences (
  run_id TEXT,
  seq_hash TEXT,
  seq_len INTEGER,
  total_count INTEGER,
  n_barcodes_present INTEGER,
  created_at TEXT,
  PRIMARY KEY (run_id, seq_hash),
  FOREIGN KEY(run_id) REFERENCES run(run_id)
);

CREATE TABLE IF NOT EXISTS sticky_sequences_by_barcode (
  run_id TEXT,
  barcode TEXT,
  seq_hash TEXT,
  count INTEGER,
  created_at TEXT,
  PRIMARY KEY (run_id, barcode, seq_hash),
  FOREIGN KEY(run_id) REFERENCES run(run_id)
);

CREATE TABLE IF NOT EXISTS sequence_clusters (
  run_id TEXT,
  seq_hash TEXT,
  cluster_id TEXT,
  seq_len INTEGER,
  created_at TEXT,
  PRIMARY KEY (run_id, seq_hash),
  FOREIGN KEY(run_id) REFERENCES run(run_id)
);

CREATE TABLE IF NOT EXISTS cluster_counts (
  run_id TEXT,
  barcode TEXT,
  cluster_id TEXT,
  count INTEGER,
  created_at TEXT,
  PRIMARY KEY (run_id, barcode, cluster_id),
  FOREIGN KEY(run_id) REFERENCES run(run_id)
);
"""


def connect_db(db_path: Path) -> sqlite3.Connection:
    ensure_parent(db_path)
    conn = sqlite3.connect(
        str(db_path),
        timeout=DB_TIMEOUT_SECONDS,
        check_same_thread=False
    )

    conn.execute("PRAGMA foreign_keys=ON;")
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute(f"PRAGMA busy_timeout={DB_BUSY_TIMEOUT_MS};")

    conn.executescript(SCHEMA_SQL)

    # Migrations for existing databases
    for ddl in [
        "CREATE TABLE IF NOT EXISTS sequence_clusters (run_id TEXT, seq_hash TEXT, cluster_id TEXT, seq_len INTEGER, created_at TEXT, PRIMARY KEY (run_id, seq_hash), FOREIGN KEY(run_id) REFERENCES run(run_id));",
        "CREATE TABLE IF NOT EXISTS cluster_counts (run_id TEXT, barcode TEXT, cluster_id TEXT, count INTEGER, created_at TEXT, PRIMARY KEY (run_id, barcode, cluster_id), FOREIGN KEY(run_id) REFERENCES run(run_id));",
    ]:
        try:
            conn.execute(ddl)
            conn.commit()
        except sqlite3.OperationalError:
            pass

    return conn


def run_exists(conn: sqlite3.Connection, run_id: str) -> bool:
    cur = conn.execute("SELECT 1 FROM run WHERE run_id=? LIMIT 1;", (run_id,))
    return cur.fetchone() is not None


def upsert_df(conn: sqlite3.Connection, df: pd.DataFrame, table: str, pk_cols: List[str]) -> None:
    """
    Upsert DataFrame rows by primary key columns.
    Works for SQLite >= 3.24 (UPSERT).
    """
    if df.empty:
        return

    cols = list(df.columns)
    placeholders = ", ".join(["?"] * len(cols))
    col_list = ", ".join(cols)

    update_cols = [c for c in cols if c not in pk_cols]
    set_clause = ", ".join([f"{c}=excluded.{c}" for c in update_cols]) if update_cols else ""

    pk_list = ", ".join(pk_cols)

    sql = f"""
    INSERT INTO {table} ({col_list})
    VALUES ({placeholders})
    ON CONFLICT ({pk_list}) DO UPDATE SET
    {set_clause};
    """.strip()

    conn.executemany(sql, df.itertuples(index=False, name=None))
    conn.commit()


def delete_run(db_path, run_id: str):
    """
    SAFE but FAST deletion - 3-5x faster than original without corruption risk.
    
    Key optimizations:
    - BEGIN IMMEDIATE (prevents lock upgrades)
    - Larger transactions (commit less often)
    - Optimized but SAFE pragmas
    - No parallel threading (avoids corruption)
    """
    # Tables in FK order (children first, parents last)
    tables = [
        "sticky_sequences_by_barcode",
        "alignment_groups_by_barcode",
        "cluster_counts",
        "sequence_clusters",
        "sticky_sequences",
        "alignment_groups",
        "barcode_qc",
        "run_qc",
        "run",
    ]

    def table_exists(conn, t):
        cur = conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?;",
            (t,),
        )
        return cur.fetchone() is not None

    def get_row_count(conn, table: str, run_id: str) -> int:
        """Fast count of rows to delete"""
        cur = conn.execute(f"SELECT COUNT(*) FROM {table} WHERE run_id=?;", (run_id,))
        return cur.fetchone()[0]

    # UI setup
    status = st.empty()
    prog = st.progress(0)
    
    # Open ONE connection for all operations
    conn = sqlite3.connect(str(db_path), timeout=60)
    
    try:
        # SAFE pragmas that improve performance without corruption risk
        conn.execute("PRAGMA busy_timeout=60000;")
        conn.execute("PRAGMA journal_mode=WAL;")  # Keep WAL (safe)
        conn.execute("PRAGMA synchronous=NORMAL;")  # NOT OFF (prevents corruption)
        conn.execute("PRAGMA cache_size=100000;")  # 100MB cache
        conn.execute("PRAGMA temp_store=MEMORY;")
        
        # Force checkpoint to see latest data
        conn.execute("PRAGMA wal_checkpoint(PASSIVE);")
        
        # Count what we're deleting
        status.write("Counting rows to delete...")
        total_rows = 0
        table_counts = {}
        
        for table in tables:
            if not table_exists(conn, table):
                continue
            count = get_row_count(conn, table, run_id)
            table_counts[table] = count
            total_rows += count
            if count > 0:
                status.write(f"  {table}: {count:,} rows")
        
        if total_rows == 0:
            status.warning(f"No data found for run_id={run_id}")
            prog.progress(100)
            return
        
        status.write(f"\nDeleting {total_rows:,} rows total...")
        deleted_so_far = 0
        
        # Delete from each table with optimized chunking
        for i, table in enumerate(tables, start=1):
            if table not in table_counts or table_counts[table] == 0:
                continue
            
            count = table_counts[table]
            status.write(f"[{i}/{len(tables)}] Deleting from {table} ({count:,} rows)...")
            
            # KEY OPTIMIZATION: Use BEGIN IMMEDIATE to avoid lock upgrade issues
            conn.execute("BEGIN IMMEDIATE;")
            
            try:
                # For small tables, delete all at once
                if count < 10_000:
                    conn.execute(f"DELETE FROM {table} WHERE run_id=?;", (run_id,))
                    deleted_so_far += count
                    conn.commit()
                else:
                    # For large tables: bigger chunks = fewer transactions = faster
                    # 100K is safe and much faster than 50K
                    chunk_size = 100_000
                    chunks_deleted = 0
                    
                    while True:
                        # Fetch rowids
                        cur = conn.execute(
                            f"SELECT rowid FROM {table} WHERE run_id=? LIMIT ?;",
                            (run_id, chunk_size),
                        )
                        rowids = [r[0] for r in cur.fetchall()]
                        
                        if not rowids:
                            break
                        
                        # Delete this chunk
                        placeholders = ",".join(["?"] * len(rowids))
                        conn.execute(f"DELETE FROM {table} WHERE rowid IN ({placeholders});", rowids)
                        
                        chunks_deleted += 1
                        deleted_so_far += len(rowids)
                        
                        # Commit every 3 chunks (300K rows) to avoid huge transactions
                        if chunks_deleted % 3 == 0:
                            conn.commit()
                            # Update progress
                            progress = int(deleted_so_far / total_rows * 100)
                            prog.progress(min(progress, 99))
                            # Start new transaction immediately
                            conn.execute("BEGIN IMMEDIATE;")
                    
                    # Final commit for this table
                    conn.commit()
                
            except Exception as e:
                conn.rollback()
                raise e
            
            progress = int(deleted_so_far / total_rows * 100)
            prog.progress(min(progress, 99))
        
        prog.progress(100)
        
        # Verify deletion
        status.write("Verifying deletion...")
        remaining = get_row_count(conn, "run", run_id)
        
        if remaining > 0:
            status.error(f"⚠️ FAILED: {remaining} rows still in 'run' table!")
        else:
            status.success(f"✅ Successfully deleted {total_rows:,} rows")
        
        # Checkpoint WAL to merge changes (makes them visible to other connections)
        status.write("Finalizing...")
        conn.execute("PRAGMA wal_checkpoint(RESTART);")
        
    except sqlite3.DatabaseError as e:
        if "malformed" in str(e):
            status.error(f"❌ DATABASE CORRUPTED: {e}")
            status.error("Your database is corrupted. You need to:")
            status.code("cd /Users/ishanghosh/Downloads/bam_pass/\nsqlite3 run_all_v2.db .dump > dump.sql\nsqlite3 run_all_v2_fixed.db < dump.sql")
        else:
            status.error(f"❌ Database error: {e}")
        raise
    finally:
        conn.close()
    
    status.success(f"✅ Complete!")


# -----------------------------
# BAM parsing - OPTIMIZED
# -----------------------------
def parse_run_meta_from_header(header_text: str) -> Dict[str, str]:
    """
    Robust ONT-ish header parsing. Never returns None; missing values are "NA".
    Uses first @RG line as canonical.
    """
    rg_lines = [ln for ln in header_text.splitlines() if ln.startswith("@RG")]
    pg_lines = [ln for ln in header_text.splitlines() if ln.startswith("@PG")]

    meta = {
        "run_id": "NA",
        "flowcell_id": "NA",
        "device_id": "NA",
        "device_type": "ONT",
        "exp_start_time": "NA",
        "sequencing_kit": "NA",
        "flowcell_type": "NA",
        "experiment_type": "NA",
        "sample_id": "NA",
    }

    def tag_value(line: str, tag: str) -> str:
        m = re.search(rf"(?:^|\t){re.escape(tag)}:([^\t]+)", line)
        return m.group(1) if m else "NA"

    if rg_lines:
        rg = rg_lines[0]
        meta["exp_start_time"] = na(tag_value(rg, "DT"))
        meta["flowcell_id"] = na(tag_value(rg, "PU"))
        meta["device_id"] = na(tag_value(rg, "PM"))
        meta["sample_id"] = na(tag_value(rg, "LB"))

        ds = tag_value(rg, "DS")
        if ds != "NA":
            m = re.search(r"runid=([a-f0-9-]+)", ds)
            if m:
                meta["run_id"] = m.group(1)
            m = re.search(r"basecall_model=([^\s]+)", ds)
            if m:
                meta["sequencing_kit"] = m.group(1)

    if pg_lines and meta["experiment_type"] == "NA":
        joined = " ".join(pg_lines)
        m = re.search(r"(?:^|\s)PN:([^\s]+)", joined)
        if m:
            meta["experiment_type"] = m.group(1)

    return meta


@dataclass
class BamFilter:
    """Configuration for BAM record filtering"""
    mapped_only: bool = False
    primary_only: bool = False
    min_mapq: int = 0
    
    def build_samtools_args(self) -> List[str]:
        """Build samtools view filter arguments"""
        args = []
        flags_exclude = 0
        
        if self.mapped_only:
            flags_exclude |= FLAG_UNMAPPED
        if self.primary_only:
            flags_exclude |= (FLAG_SECONDARY | FLAG_SUPPLEMENTARY)
            
        if flags_exclude:
            args += ["-F", str(flags_exclude)]
        if self.min_mapq > 0:
            args += ["-q", str(self.min_mapq)]
            
        return args


@dataclass
class BamRecord:
    """Minimal parsed BAM record for stats collection"""
    qname: str
    flag: int
    rname: str
    mapq: int
    seq: str
    
    @property
    def is_mapped(self) -> bool:
        return not (self.flag & FLAG_UNMAPPED)
    
    @property
    def is_primary(self) -> bool:
        return not (self.flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY))
    
    @property
    def is_secondary(self) -> bool:
        return bool(self.flag & FLAG_SECONDARY)
    
    @property
    def is_supplementary(self) -> bool:
        return bool(self.flag & FLAG_SUPPLEMENTARY)
    
    @property
    def seq_len(self) -> int:
        return 0 if self.seq == "*" else len(self.seq)


def stream_bam_records(bam_path: Path, filters: Optional[BamFilter] = None) -> Iterator[BamRecord]:
    """
    Unified BAM streaming with configurable filters.
    Single source of truth for BAM iteration.
    """
    view_cmd = ["samtools", "view"]
    
    if filters:
        view_cmd += filters.build_samtools_args()
    
    view_cmd += [str(bam_path)]
    
    p = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    assert p.stdout is not None
    for line in p.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 11:
            continue
        
        try:
            yield BamRecord(
                qname=parts[0],
                flag=int(parts[1]),
                rname=parts[2],
                mapq=int(parts[4]),
                seq=parts[9],
            )
        except (ValueError, IndexError):
            continue
    
    stderr = p.stderr.read() if p.stderr else ""
    rc = p.wait()
    if rc != 0:
        raise RuntimeError(f"samtools view failed for {bam_path}\n{stderr}")


@dataclass
class BamStats:
    """Complete BAM statistics"""
    total_reads: int
    mapped_reads: int
    mean_read_len: float
    mapq_mean: float
    n_primary: int
    n_secondary: int
    n_supplementary: int
    
    # Per-reference alignment stats
    ref_counts: Dict[str, int]
    ref_mapq_sums: Dict[str, float]
    ref_mapq_counts: Dict[str, int]
    
    def ref_stats(self) -> List[Tuple[str, int, float]]:
        """Get (ref_name, read_count, mapq_mean) sorted by read count"""
        result = []
        for ref, count in self.ref_counts.items():
            mean_mq = (self.ref_mapq_sums[ref] / self.ref_mapq_counts[ref]) if self.ref_mapq_counts[ref] > 0 else 0.0
            result.append((ref, count, float(mean_mq)))
        result.sort(key=lambda x: x[1], reverse=True)
        return result


def compute_bam_stats_single_pass(bam_path: Path) -> BamStats:
    """
    OPTIMIZED: Single-pass BAM statistics collection.
    Collects all stats (counts, lengths, MAPQ, per-ref) in one streaming pass.
    This replaces multiple separate samtools calls.
    """
    # Stream all records once with no filters
    total_reads = 0
    mapped_reads = 0
    n_primary = 0
    n_secondary = 0
    n_supplementary = 0
    
    sum_read_len = 0.0
    sum_mapq = 0.0
    n_for_means = 0
    
    ref_counts: Dict[str, int] = defaultdict(int)
    ref_mapq_sums: Dict[str, float] = defaultdict(float)
    ref_mapq_counts: Dict[str, int] = defaultdict(int)
    
    records_seen = 0
    for record in stream_bam_records(bam_path, filters=None):
        records_seen += 1
        total_reads += 1
        
        if record.is_mapped:
            mapped_reads += 1
        
        if record.is_primary:
            n_primary += 1
            sum_read_len += record.seq_len
            sum_mapq += record.mapq
            n_for_means += 1
            
            # Per-reference stats (only for primary mapped reads)
            if record.is_mapped and record.rname != "*":
                ref_counts[record.rname] += 1
                ref_mapq_sums[record.rname] += record.mapq
                ref_mapq_counts[record.rname] += 1
        
        elif record.is_secondary:
            n_secondary += 1
        elif record.is_supplementary:
            n_supplementary += 1
    
    # If no records seen, BAM is likely empty
    if records_seen == 0:
        pass  # stats will all be 0, which is correct

    mean_len = (sum_read_len / n_for_means) if n_for_means > 0 else 0.0
    mean_mapq = (sum_mapq / n_for_means) if n_for_means > 0 else 0.0
    
    return BamStats(
        total_reads=total_reads,
        mapped_reads=mapped_reads,
        mean_read_len=float(mean_len),
        mapq_mean=float(mean_mapq),
        n_primary=n_primary,
        n_secondary=n_secondary,
        n_supplementary=n_supplementary,
        ref_counts=dict(ref_counts),
        ref_mapq_sums=dict(ref_mapq_sums),
        ref_mapq_counts=dict(ref_mapq_counts),
    )


# -----------------------------
# Sticky sequences - OPTIMIZED
# -----------------------------
def seq_hash_sha256(seq: str) -> str:
    """Use SHA256 instead of MD5 for better collision resistance"""
    return hashlib.sha256(seq.encode("utf-8")).hexdigest()[:32]  # truncate for storage


def build_sticky_sequences_batch(
    run_id: str,
    bam_paths: List[Tuple[Path, str]],  # [(bam_path, barcode), ...]
    conn: sqlite3.Connection,
    min_mapq: int = 10,
    max_unique: Optional[int] = None,
    ui_status=None,
    ui_progress=None,
) -> None:
    """
    OPTIMIZED: Build sticky sequences with batched DB inserts.
    Processes all BAMs once and inserts in batches to reduce DB overhead.
    
    Args:
        bam_paths: List of (bam_path, barcode) tuples
    """
    created_at = utc_iso()
    
    # Accumulate all data in memory first
    per_bc: Dict[Tuple[str, str], int] = defaultdict(int)  # (barcode, hash) -> count
    overall: Dict[str, Tuple[int, int]] = {}  # hash -> (seq_len, total_count)
    
    filters = BamFilter(mapped_only=True, primary_only=True, min_mapq=min_mapq)
    
    total_bams = len(bam_paths)
    for i, (bam_path, barcode) in enumerate(bam_paths, start=1):
        if ui_status is not None:
            ui_status.write(f"Sticky scan {i}/{total_bams}: {barcode}")
        if ui_progress is not None:
            ui_progress.progress(int((i - 1) / max(1, total_bams) * 100))
        
        for record in stream_bam_records(bam_path, filters):
            if not record.seq or record.seq == "*":
                continue
                
            h = seq_hash_sha256(record.seq)
            per_bc[(barcode, h)] += 1
            
            if h in overall:
                sl, tot = overall[h]
                overall[h] = (sl, tot + 1)
            else:
                overall[h] = (len(record.seq), 1)
            
            if max_unique is not None and len(overall) >= max_unique:
                break
        
        if max_unique is not None and len(overall) >= max_unique:
            break
    
    if ui_progress is not None:
        ui_progress.progress(100)
    
    # Calculate barcode presence
    bc_set_by_hash: Dict[str, set] = defaultdict(set)
    for (bc, h) in per_bc.keys():
        bc_set_by_hash[h].add(bc)
    
    # Batch insert sticky_sequences
    if ui_status:
        ui_status.write("Inserting sticky sequences...")
    
    sticky_rows = [
        (run_id, h, int(sl), int(tot), int(len(bc_set_by_hash[h])), created_at)
        for h, (sl, tot) in overall.items()
    ]
    
    # Insert in batches to avoid memory issues with huge datasets
    for i in range(0, len(sticky_rows), STICKY_BATCH_SIZE):
        batch = sticky_rows[i:i + STICKY_BATCH_SIZE]
        df = pd.DataFrame(batch, columns=["run_id", "seq_hash", "seq_len", "total_count", "n_barcodes_present", "created_at"])
        upsert_df(conn, df, "sticky_sequences", pk_cols=["run_id", "seq_hash"])
    
    # Batch insert sticky_sequences_by_barcode
    if ui_status:
        ui_status.write("Inserting sticky sequences by barcode...")
    
    sticky_bc_rows = [
        (run_id, bc, h, int(c), created_at)
        for (bc, h), c in per_bc.items()
    ]
    
    for i in range(0, len(sticky_bc_rows), STICKY_BATCH_SIZE):
        batch = sticky_bc_rows[i:i + STICKY_BATCH_SIZE]
        df = pd.DataFrame(batch, columns=["run_id", "barcode", "seq_hash", "count", "created_at"])
        upsert_df(conn, df, "sticky_sequences_by_barcode", pk_cols=["run_id", "barcode", "seq_hash"])


# -----------------------------
# Protein Translation
# -----------------------------
def translate_dna_to_protein(dna_seq: str, frame: int = 0) -> str:
    """
    Translate DNA to protein in specified frame.
    Returns empty string if contains stop codon or invalid.
    """
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    
    dna_seq = dna_seq.upper()[frame:]
    protein = []
    
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        if len(codon) != 3:
            break
        aa = codon_table.get(codon, 'X')
        if aa == '*':  # stop codon
            break
        protein.append(aa)
    
    return ''.join(protein)


def find_best_orf(dna_seq: str, min_aa_len: int = 100) -> Optional[str]:
    """
    Find the longest ORF in all 3 forward frames.
    Returns protein sequence or None if no valid ORF found.
    Nanobodies are ~110-140 AA, so min 100 AA is reasonable.
    """
    best_protein = None
    best_len = 0
    
    for frame in range(3):
        protein = translate_dna_to_protein(dna_seq, frame)
        if len(protein) >= min_aa_len and len(protein) > best_len:
            best_protein = protein
            best_len = len(protein)
    
    return best_protein


# -----------------------------
# MMseqs2 Sequence Clustering
# -----------------------------
MIN_SEQ_LEN = 400
MAX_SEQ_LEN = 700
MMSEQS_MIN_SEQ_ID = 0.95
MMSEQS_COVERAGE = 0.9


def check_mmseqs() -> None:
    try:
        result = subprocess.run(["mmseqs", "version"], capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError("mmseqs version check failed")
    except FileNotFoundError:
        st.error("MMseqs2 is required but not found on PATH. Install with: conda install -c conda-forge -c bioconda mmseqs2")
        st.stop()


def extract_sequences_from_bam(
    bam_path: Path,
    barcode: str,
    min_len: int = MIN_SEQ_LEN,
    max_len: int = MAX_SEQ_LEN,
    min_mapq: int = 10,
) -> Dict[str, str]:
    """
    Extract filtered sequences from BAM as {seq_hash: sequence}.
    Filters: primary, mapped, MAPQ >= min_mapq, length in [min_len, max_len].
    """
    filters = BamFilter(mapped_only=True, primary_only=True, min_mapq=min_mapq)
    seqs = {}
    for record in stream_bam_records(bam_path, filters):
        if record.seq == "*" or not record.seq:
            continue
        if not (min_len <= len(record.seq) <= max_len):
            continue
        h = seq_hash_sha256(record.seq)
        seqs[h] = record.seq
    return seqs


def run_mmseqs_cluster(
    sequences: Dict[str, str],
    work_dir: Path,
    seq_counts: Optional[Dict[str, int]] = None,
    min_seq_id: float = MMSEQS_MIN_SEQ_ID,
    coverage: float = MMSEQS_COVERAGE,
    max_seqs: Optional[int] = None,
    is_protein: bool = True,
) -> Dict[str, str]:
    """
    Cluster sequences with MMseqs2 easy-linclust.
    Returns {seq_hash: cluster_rep_hash}.
    
    max_seqs: Max sequences to cluster. If None, auto-determined by available RAM.
    is_protein: If True, uses 80% identity (proteins cluster better).
                If False, uses min_seq_id parameter (nucleotides).
    """
    if max_seqs is None:
        max_seqs = recommend_max_sequences()
    if not sequences:
        return {}

    cluster_dir = work_dir / "mmseqs_cluster"
    cluster_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = cluster_dir / "input.fasta"
    result_prefix = str(cluster_dir / "clusters")
    tmp_dir = cluster_dir / "tmp"
    tmp_dir.mkdir(exist_ok=True)

    # Sort by count descending so we keep the most abundant sequences
    if seq_counts:
        sorted_seqs = sorted(sequences.items(), key=lambda x: seq_counts.get(x[0], 1), reverse=True)
    else:
        sorted_seqs = list(sequences.items())

    if len(sorted_seqs) > max_seqs:
        sorted_seqs = sorted_seqs[:max_seqs]

    with open(fasta_path, "w") as f:
        for h, seq in sorted_seqs:
            f.write(f">{h}\n{seq}\n")

    # Proteins cluster much better at 80% identity (accounts for sequencing errors + real variation)
    actual_min_id = 0.80 if is_protein else min_seq_id
    
    run_cmd([
        "mmseqs", "easy-linclust",
        str(fasta_path),
        result_prefix,
        str(tmp_dir),
        "--min-seq-id", str(actual_min_id),
        "-c", str(coverage),
        "--cov-mode", "0",
        "--threads", str(MAX_WORKERS),
        "--split-memory-limit", "3G",
    ])

    # Parse cluster TSV: rep_hash \t member_hash
    cluster_tsv = Path(result_prefix + "_cluster.tsv")
    hash_to_cluster: Dict[str, str] = {}

    if cluster_tsv.exists():
        with open(cluster_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    rep_hash, member_hash = parts
                    hash_to_cluster[member_hash] = rep_hash

    # Cleanup tmp
    shutil.rmtree(str(tmp_dir), ignore_errors=True)

    return hash_to_cluster


def build_sequence_clusters(
    run_id: str,
    bam_paths: List[Tuple[Path, str]],
    conn: sqlite3.Connection,
    work_dir: Path,
    min_mapq: int = 10,
    min_len: int = None,
    max_len: int = None,
    max_seqs: Optional[int] = None,
    ui_status=None,
    ui_progress=None,
) -> None:
    """
    Full pipeline: extract sequences → translate to protein → cluster with MMseqs2 → store in DB.
    Clusters across ALL barcodes combined, then counts per barcode.
    
    min_len, max_len: If None, uses MIN_SEQ_LEN, MAX_SEQ_LEN constants
    max_seqs: If None, auto-determined by RAM
    """
    if min_len is None:
        min_len = MIN_SEQ_LEN
    if max_len is None:
        max_len = MAX_SEQ_LEN
    
    created_at = utc_iso()
    total_bams = len(bam_paths)

    # Step 1: Extract sequences from all barcodes
    if ui_status:
        ui_status.write("Extracting sequences from BAMs...")

    all_seqs: Dict[str, str] = {}  # hash -> seq (deduplicated across barcodes)
    per_barcode_counts: Dict[str, Dict[str, int]] = {}  # barcode -> {hash: count}
    global_counts: Dict[str, int] = defaultdict(int)  # hash -> total count across all barcodes

    for i, (bam_path, barcode) in enumerate(bam_paths, start=1):
        if ui_status:
            ui_status.write(f"Extracting {i}/{total_bams}: {barcode}")
        if ui_progress:
            ui_progress.progress(int((i - 1) / total_bams * 50))

        # Single pass: collect DNAseqs, translate to protein, count per barcode
        filters = BamFilter(mapped_only=True, primary_only=True, min_mapq=min_mapq)
        bc_counts: Dict[str, int] = defaultdict(int)
        for record in stream_bam_records(bam_path, filters):
            if not record.seq or record.seq == "*":
                continue
            if not (min_len <= len(record.seq) <= max_len):
                continue
            
            # Translate DNA to protein
            protein = find_best_orf(record.seq, min_aa_len=100)
            if not protein:
                continue  # Skip if no valid ORF
            
            # Hash the protein sequence
            h = seq_hash_sha256(protein)
            all_seqs[h] = protein
            bc_counts[h] += 1
            global_counts[h] += 1

        per_barcode_counts[barcode] = dict(bc_counts)

    if not all_seqs:
        if ui_status:
            ui_status.error("No sequences extracted. Check length filter and MAPQ settings.")
        return

    # Determine actual cap
    actual_cap = max_seqs if max_seqs is not None else recommend_max_sequences()
    
    if ui_status:
        ui_status.write(f"Extracted {len(all_seqs):,} unique protein sequences across all barcodes. Clustering top {actual_cap:,} at 80% identity...")
    if ui_progress:
        ui_progress.progress(50)

    # Step 2: Cluster proteins at 80% identity (much more robust than DNA at 95%)
    hash_to_cluster = run_mmseqs_cluster(all_seqs, work_dir, seq_counts=dict(global_counts), max_seqs=max_seqs, is_protein=True)

    if ui_status:
        n_clusters = len(set(hash_to_cluster.values()))
        ui_status.write(f"MMseqs2 found {n_clusters:,} clusters from top {actual_cap:,} sequences (out of {len(all_seqs):,} total).")
    if ui_progress:
        ui_progress.progress(75)

    # Step 3: Store sequence_clusters rows
    if ui_status:
        ui_status.write("Storing cluster assignments...")

    sc_rows = []
    for h, seq in all_seqs.items():
        cluster_id = hash_to_cluster.get(h, h)  # fallback to self if not in output
        sc_rows.append((run_id, h, cluster_id, len(seq), created_at))

    # Insert in batches - note: no barcode column here, it's per-sequence
    for i in range(0, len(sc_rows), STICKY_BATCH_SIZE):
        batch = sc_rows[i:i + STICKY_BATCH_SIZE]
        df = pd.DataFrame(batch, columns=["run_id", "seq_hash", "cluster_id", "seq_len", "created_at"])
        upsert_df(conn, df, "sequence_clusters", pk_cols=["run_id", "seq_hash"])

    # Step 4: Build cluster_counts per barcode
    if ui_status:
        ui_status.write("Computing cluster counts per barcode...")

    cc_rows = []
    for barcode, bc_counts in per_barcode_counts.items():
        cluster_counter: Dict[str, int] = defaultdict(int)
        for h, count in bc_counts.items():
            cluster_id = hash_to_cluster.get(h, h)
            cluster_counter[cluster_id] += count
        for cluster_id, count in cluster_counter.items():
            cc_rows.append((run_id, barcode, cluster_id, count, created_at))

    for i in range(0, len(cc_rows), STICKY_BATCH_SIZE):
        batch = cc_rows[i:i + STICKY_BATCH_SIZE]
        df = pd.DataFrame(batch, columns=["run_id", "barcode", "cluster_id", "count", "created_at"])
        upsert_df(conn, df, "cluster_counts", pk_cols=["run_id", "barcode", "cluster_id"])

    if ui_progress:
        ui_progress.progress(100)
    if ui_status:
        ui_status.success(f"✅ Clustering complete! {len(set(hash_to_cluster.values())):,} clusters stored.")


def generate_sql_from_natural_language(question: str, schema: str) -> str:
    """
    Use Claude API to convert natural language question to SQL query.
    Returns the SQL query as a string.
    """
    client = anthropic.Anthropic(api_key=ANTHROPIC_API_KEY)
    
    prompt = f"""You are a SQL expert helping biologists query their MinION nanobody sequencing database.

DATABASE SCHEMA:
{schema}

USER QUESTION: {question}

Generate a SQL query to answer this question. Return ONLY the SQL query, no explanation.
The query should be safe (SELECT only, no DELETE/UPDATE/DROP).
If the question is unclear, make reasonable assumptions about what the user wants.

SQL Query:"""

    message = client.messages.create(
        model="claude-3-haiku-20240307",
        max_tokens=1024,
        messages=[{"role": "user", "content": prompt}]
    )
    
    sql = message.content[0].text.strip()
    
    # Remove markdown code fences if present
    sql = sql.replace("```sql", "").replace("```", "").strip()
    
    return sql


# -----------------------------
# Zip ingest - OPTIMIZED
# -----------------------------
@dataclass
class ZipBamEntry:
    bam_member: zipfile.ZipInfo
    bai_member: Optional[zipfile.ZipInfo]
    barcode: str


def find_bams_in_zip(zippath: Path) -> List[ZipBamEntry]:
    with zipfile.ZipFile(zippath, "r") as z:
        members = z.infolist()

        bam_infos = [m for m in members if m.filename.lower().endswith(".bam") and not m.filename.lower().endswith(".bam.bai")]
        bai_by_base: Dict[str, zipfile.ZipInfo] = {}
        for m in members:
            if m.filename.lower().endswith(".bam.bai"):
                base = m.filename[:-4]
                bai_by_base[base] = m

        out: List[ZipBamEntry] = []
        for bam in bam_infos:
            base = bam.filename
            bai = bai_by_base.get(base, None)

            parts = Path(bam.filename).parts
            barcode = "unknown"
            if len(parts) >= 2:
                barcode = parts[-2]
            out.append(ZipBamEntry(bam_member=bam, bai_member=bai, barcode=barcode))

        out.sort(key=lambda e: (e.barcode, e.bam_member.filename))
        return out


def safe_extract_member(z: zipfile.ZipFile, member: zipfile.ZipInfo, dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    out_path = dest / Path(member.filename).name
    with z.open(member, "r") as src, open(out_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out_path


def process_single_bam_simple(
    zip_path_str: str,
    bam_filename: str,
    bai_filename: Optional[str],
    barcode: str,
    work_dir_str: str,
    run_id: str,
    created_at: str,
) -> Dict:
    """
    Process a single BAM file (for parallel execution).
    Uses only simple types that can be pickled.
    """
    try:
        import zipfile
        from pathlib import Path
        import shutil
        import os
        from collections import defaultdict
        
        zip_path = Path(zip_path_str)
        work_dir = Path(work_dir_str)
        
        # Extract BAM to work_dir
        with zipfile.ZipFile(zip_path, "r") as z:
            tmp_dir = work_dir / f"tmp_{barcode}_{os.getpid()}"
            tmp_dir.mkdir(parents=True, exist_ok=True)
            
            try:
                # Extract BAM
                bam_member = z.getinfo(bam_filename)
                bam_dest = tmp_dir / Path(bam_filename).name
                with z.open(bam_member, "r") as src, open(bam_dest, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                
                # Extract BAI if exists
                if bai_filename:
                    bai_member = z.getinfo(bai_filename)
                    bai_dest = tmp_dir / Path(bai_filename).name
                    with z.open(bai_member, "r") as src, open(bai_dest, "wb") as dst:
                        shutil.copyfileobj(src, dst)
                
                # Get header for metadata
                header = run_cmd(["samtools", "view", "-H", str(bam_dest)]).stdout
                meta = parse_run_meta_from_header(header)
                
                # Single-pass stats collection
                stats = compute_bam_stats_single_pass(bam_dest)
                
                mapping_pct = (100.0 * stats.mapped_reads / stats.total_reads) if stats.total_reads > 0 else 0.0
                
                barcode_row = (
                    run_id,
                    barcode,
                    int(stats.total_reads),
                    int(stats.mapped_reads),
                    float(mapping_pct),
                    float(stats.mean_read_len),
                    float(stats.mapq_mean),
                    created_at,
                )
                
                # Per-reference stats
                agb_rows = [
                    (run_id, barcode, ref_name, int(read_count), float(mapq_mean), created_at)
                    for ref_name, read_count, mapq_mean in stats.ref_stats()
                ]
                
                # Convert stats to dict for pickling
                stats_dict = {
                    'total_reads': stats.total_reads,
                    'mapped_reads': stats.mapped_reads,
                    'mean_read_len': stats.mean_read_len,
                    'mapq_mean': stats.mapq_mean,
                    'n_primary': stats.n_primary,
                    'n_secondary': stats.n_secondary,
                    'n_supplementary': stats.n_supplementary,
                    'ref_counts': stats.ref_counts,
                    'ref_mapq_sums': stats.ref_mapq_sums,
                    'ref_mapq_counts': stats.ref_mapq_counts,
                }
                
                return {
                    'barcode': barcode,
                    'barcode_row': barcode_row,
                    'agb_rows': agb_rows,
                    'stats': stats_dict,
                    'meta': meta,
                }
            
            finally:
                # Cleanup
                if tmp_dir.exists():
                    shutil.rmtree(tmp_dir, ignore_errors=True)
    
    except Exception as e:
        import traceback
        return {
            'barcode': barcode,
            'error': str(e),
            'traceback': traceback.format_exc(),
        }


def ingest_zip_into_db(
    zippath: Path,
    conn: sqlite3.Connection,
    work_dir: Path,
    skip_if_run_exists: bool,
    ui_status,
    ui_progress,
) -> Optional[str]:
    """
    Ingest zip with parallel BAM processing.
    """
    check_samtools()
    work_dir.mkdir(parents=True, exist_ok=True)

    ui_status.write(f"Zip: {zippath}")
    entries = find_bams_in_zip(zippath)
    if not entries:
        ui_status.error("No BAMs found in zip.")
        return None

    ui_status.write(f"Found {len(entries)} BAM(s)")

    # Extract first BAM to get run_id
    with zipfile.ZipFile(zippath, "r") as z:
        first = entries[0]
        tmp_dir = work_dir / "tmp_first"
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir, ignore_errors=True)
        tmp_dir.mkdir(parents=True, exist_ok=True)

        bam_path = safe_extract_member(z, first.bam_member, tmp_dir)
        if first.bai_member is not None:
            safe_extract_member(z, first.bai_member, tmp_dir)

        header = run_cmd(["samtools", "view", "-H", str(bam_path)]).stdout
        meta = parse_run_meta_from_header(header)
        run_id = meta["run_id"]

        shutil.rmtree(tmp_dir, ignore_errors=True)

    if run_id == "NA":
        run_id = f"run_{hashlib.md5(str(zippath).encode()).hexdigest()[:10]}"
        ui_status.warning(f"Could not parse run_id from BAM header; using {run_id}")

    if skip_if_run_exists and run_exists(conn, run_id):
        ui_status.warning(f"Run already exists in DB: {run_id} (skipping ingest).")
        return None

    created_at = utc_iso()

    # Aggregate statistics
    total_reads = 0
    mapped_reads = 0
    sum_read_len = 0.0
    sum_mapq = 0.0
    n_for_means = 0

    n_primary = 0
    n_secondary = 0
    n_supp = 0

    barcode_rows = []
    agb_rows = []

    ref_count_total: Dict[str, int] = defaultdict(int)
    ref_mapq_sum: Dict[str, float] = defaultdict(float)
    ref_mapq_n: Dict[str, int] = defaultdict(int)

    total_items = len(entries)
    ui_progress.progress(0)

    ui_status.write(f"Processing {len(entries)} BAMs in parallel (max {MAX_WORKERS} workers)...")
    
    results = []
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(
                process_single_bam_simple,
                str(zippath),
                entry.bam_member.filename,
                entry.bai_member.filename if entry.bai_member else None,
                entry.barcode,
                str(work_dir),
                run_id,
                created_at
            ): i
            for i, entry in enumerate(entries, start=1)
        }
        
        for future in as_completed(futures):
            i = futures[future]
            try:
                result = future.result()
                
                if 'error' in result:
                    ui_status.error(f"❌ BAM {i} ({result['barcode']}) failed: {result['error']}")
                    continue
                
                results.append(result)
                ui_status.write(f"Processed {len(results)}/{total_items}: {result['barcode']}")
                ui_progress.progress(int(len(results) / max(1, total_items) * 100))
                
            except Exception as e:
                ui_status.error(f"Error processing BAM {i}: {str(e)}")
    
    for result in results:
        stats_dict = result['stats']
        
        total_reads += stats_dict['total_reads']
        mapped_reads += stats_dict['mapped_reads']
        sum_read_len += stats_dict['mean_read_len'] * max(1, stats_dict['n_primary'])
        sum_mapq += stats_dict['mapq_mean'] * max(1, stats_dict['n_primary'])
        n_for_means += max(1, stats_dict['n_primary'])
        
        n_primary += stats_dict['n_primary']
        n_secondary += stats_dict['n_secondary']
        n_supp += stats_dict['n_supplementary']
        
        barcode_rows.append(result['barcode_row'])
        agb_rows.extend(result['agb_rows'])
        
        for k in ["flowcell_id", "device_id", "exp_start_time", "sequencing_kit", "sample_id", "experiment_type"]:
            if meta.get(k, "NA") == "NA" and result['meta'].get(k, "NA") != "NA":
                meta[k] = result['meta'][k]
        
        for ref_name, count in stats_dict['ref_counts'].items():
            ref_count_total[ref_name] += count
            ref_mapq_sum[ref_name] += (stats_dict['ref_mapq_sums'][ref_name] / stats_dict['ref_mapq_counts'][ref_name]) * count if stats_dict['ref_mapq_counts'][ref_name] > 0 else 0
            ref_mapq_n[ref_name] += count

    # Build DataFrames and insert
    run_df = pd.DataFrame(
        [
            (
                run_id,
                meta["flowcell_id"],
                meta["device_id"],
                meta["device_type"],
                meta["exp_start_time"],
                meta["sequencing_kit"],
                meta["flowcell_type"],
                meta["experiment_type"],
                meta["sample_id"],
                created_at,
            )
        ],
        columns=[
            "run_id",
            "flowcell_id",
            "device_id",
            "device_type",
            "exp_start_time",
            "sequencing_kit",
            "flowcell_type",
            "experiment_type",
            "sample_id",
            "created_at",
        ],
    )

    mean_read_len = (sum_read_len / n_for_means) if n_for_means else 0.0
    mapping_pct_total = (100.0 * mapped_reads / total_reads) if total_reads else 0.0

    run_qc_df = pd.DataFrame(
        [
            (
                run_id,
                int(total_reads),
                float(mean_read_len),
                int(n_primary),
                int(n_secondary),
                int(n_supp),
                int(mapped_reads),
                float(mapping_pct_total),
                created_at,
            )
        ],
        columns=[
            "run_id",
            "n_reads_total",
            "mean_read_length",
            "n_primary",
            "n_secondary",
            "n_supplementary",
            "mapped_reads_total",
            "mapping_pct",
            "created_at",
        ],
    )

    barcode_qc_df = pd.DataFrame(
        barcode_rows,
        columns=[
            "run_id",
            "barcode",
            "total_reads",
            "mapped_reads_total",
            "mapping_pct",
            "mean_read_length",
            "mapq_mean",
            "created_at",
        ],
    )

    ag_rows = []
    for ref_name, cnt in sorted(ref_count_total.items(), key=lambda x: x[1], reverse=True):
        mean_mq = (ref_mapq_sum[ref_name] / ref_mapq_n[ref_name]) if ref_mapq_n.get(ref_name, 0) else 0.0
        ag_rows.append((run_id, ref_name, int(cnt), float(mean_mq), created_at))

    alignment_groups_df = pd.DataFrame(ag_rows, columns=["run_id", "ref_name", "read_count", "mapq_mean", "created_at"])
    agb_df = pd.DataFrame(agb_rows, columns=["run_id", "barcode", "ref_name", "read_count", "mapq_mean", "created_at"])

    upsert_df(conn, run_df, "run", pk_cols=["run_id"])
    upsert_df(conn, run_qc_df, "run_qc", pk_cols=["run_id"])
    upsert_df(conn, barcode_qc_df, "barcode_qc", pk_cols=["run_id", "barcode"])
    upsert_df(conn, alignment_groups_df, "alignment_groups", pk_cols=["run_id", "ref_name"])
    upsert_df(conn, agb_df, "alignment_groups_by_barcode", pk_cols=["run_id", "barcode", "ref_name"])

    ui_status.success(f"Ingest complete: run_id={run_id}")
    ui_progress.progress(100)
    return run_id


# -----------------------------
# UI helpers
# -----------------------------
def download_button_df(label: str, df: pd.DataFrame, filename: str):
    st.download_button(
        label=label,
        data=df.to_csv(index=False).encode("utf-8"),
        file_name=filename,
        mime="text/csv",
        use_container_width=True,
    )


def sql_df(conn: sqlite3.Connection, query: str, params: Tuple = ()) -> pd.DataFrame:
    return pd.read_sql_query(query, conn, params=params)


def run_selectbox(conn: sqlite3.Connection, label: str, key: str) -> Optional[str]:
    """Selectbox showing run_ids."""
    runs = sql_df(conn, "SELECT run_id FROM run ORDER BY exp_start_time DESC, created_at DESC;")
    if runs.empty:
        return None
    options = runs["run_id"].tolist()
    return st.selectbox(label, options, key=key)


# -----------------------------
# Pages
# -----------------------------
def page_overview(conn: sqlite3.Connection):
    st.header("Overview")

    runs = sql_df(conn, "SELECT * FROM run ORDER BY exp_start_time DESC, created_at DESC;")
    if runs.empty:
        st.info("No runs in the database yet. Go to **Ingest** to add one.")
        return

    st.subheader("Runs")
    st.dataframe(runs, use_container_width=True)
    download_button_df("Download Runs CSV", runs, "runs.csv")

    run_qc = sql_df(conn, "SELECT * FROM run_qc ORDER BY created_at DESC;")
    st.subheader("Run QC")
    st.dataframe(run_qc, use_container_width=True)
    download_button_df("Download Run QC CSV", run_qc, "run_qc.csv")


def page_run_detail(conn: sqlite3.Connection):
    st.header("Run detail")

    if sql_df(conn, "SELECT 1 FROM run LIMIT 1;").empty:
        st.info("No runs in the database yet.")
        return

    run_id = run_selectbox(conn, "Select run", key="run_detail_select")
    if run_id is None:
        return

    run_row = sql_df(conn, """
        SELECT run_id, flowcell_id, device_id, device_type, exp_start_time,
               sequencing_kit, flowcell_type, experiment_type, sample_id, created_at
        FROM run WHERE run_id=?;
    """, (run_id,))
    run_qc = sql_df(conn, "SELECT * FROM run_qc WHERE run_id=?;", (run_id,))
    st.subheader("Run metadata")
    st.dataframe(run_row, use_container_width=True)

    st.subheader("Run QC")
    st.dataframe(run_qc, use_container_width=True)

    st.divider()

    barcode_qc = sql_df(conn, "SELECT * FROM barcode_qc WHERE run_id=? ORDER BY mapping_pct DESC;", (run_id,))
    st.subheader("Barcode QC (ranked by mapping_pct)")
    st.dataframe(barcode_qc, use_container_width=True)
    download_button_df("Download Barcode QC CSV", barcode_qc, f"{run_id}_barcode_qc.csv")

    if not barcode_qc.empty:
        selected_barcode = st.selectbox("Inspect barcode", barcode_qc["barcode"].tolist())

        barcode_align = sql_df(
            conn,
            """
            SELECT ref_name, read_count, mapq_mean
            FROM alignment_groups_by_barcode
            WHERE run_id=? AND barcode=?
            ORDER BY read_count DESC;
            """,
            (run_id, selected_barcode),
        )

        st.subheader(f"Alignment for {selected_barcode}")
        st.dataframe(barcode_align, use_container_width=True)
        download_button_df("Download Alignment CSV", barcode_align, f"{run_id}_{selected_barcode}_alignment.csv")

    st.divider()

    align_total = sql_df(conn, "SELECT * FROM alignment_groups WHERE run_id=? ORDER BY read_count DESC;", (run_id,))
    st.subheader("Alignment groups (overall)")
    st.dataframe(align_total, use_container_width=True)
    download_button_df("Download Overall Alignment CSV", align_total, f"{run_id}_alignment_groups.csv")


def page_ingest(conn: sqlite3.Connection, db_path: Path):
    st.header("Ingest")

    # --- Delete run (in-app) ---
    st.subheader("Delete run (for re-ingest)")
    
    runs = sql_df(conn, "SELECT run_id, created_at FROM run ORDER BY created_at DESC;")
    if runs.empty:
        st.info("No runs in this database yet.")
    else:
        run_to_delete = run_selectbox(conn, "Select run to delete", key="delete_run_id")
        confirm = st.checkbox("I understand this permanently deletes this run from the selected DB", value=False)
        col_del1, _ = st.columns([1, 2])
        with col_del1:
            if st.button("Delete run", type="secondary", disabled=not confirm, use_container_width=True):
                try:
                    # Close the main connection to release any locks
                    conn.close()
                    time.sleep(0.1)
                    
                    # Perform deletion
                    delete_run(db_path, run_to_delete)
                    
                    # Clear cache and rerun
                    st.cache_data.clear()
                    
                except Exception as e:
                    st.error(f"Deletion failed: {e}")
                    import traceback
                    st.code(traceback.format_exc())
                    st.stop()
                
                # Rerun to refresh with new connection
                st.rerun()

    st.divider()

    zip_path_disk = st.text_input("Zip path on disk", value="", placeholder="/Users/you/Downloads/run4.zip")

    st.subheader("Ingest settings")
    work_dir = st.text_input(
        "Work dir (temporary extraction per BAM)",
        value=str(Path(tempfile.gettempdir()) / "minion_tmp"),
        help="Use an external drive if your system disk is low on space.",
    )
    work_dir = Path(work_dir)

    skip_if_exists = st.checkbox("Skip ingest if run_id already exists in DB", value=True)

    st.divider()
    st.subheader("Database")
    st.write(f"DB: `{db_path}`")

    zip_path: Optional[Path] = None

    if zip_path_disk.strip():
        zp = Path(zip_path_disk.strip())
        if not zp.exists():
            st.error("Zip path does not exist.")
        else:
            zip_path = zp

    if zip_path is None:
        st.info("Enter a zip path above to get started.")
        return

    if st.button("Ingest now", type="primary", use_container_width=True):
        # Remember for sticky sequences page
        st.session_state["last_zip_path"] = str(zip_path)
        st.session_state["last_work_dir"] = str(work_dir)

        status = st.empty()
        progress = st.progress(0)
        try:
            rid = ingest_zip_into_db(
                zippath=zip_path,
                conn=conn,
                work_dir=work_dir,
                skip_if_run_exists=skip_if_exists,
                ui_status=status,
                ui_progress=progress,
            )
            if rid:
                st.success(f"Done. Ingested run_id={rid}")
            else:
                st.warning("No run ingested (likely skipped because it already exists).")
        finally:
            pass


def page_query(conn: sqlite3.Connection):
    st.header("SQL Query")

    st.caption("Ask questions in plain English or write SQL queries directly. Results can be downloaded as CSV.")

    # AI Chatbot section
    st.subheader("🤖 Ask Claude (Natural Language)")
    st.caption("Type your question in plain English and Claude will generate the SQL query for you.")
    
    user_question = st.text_input(
        "What do you want to know?",
        placeholder="e.g., Show me the top 10 clusters in barcode08 with at least 50 reads",
        key="ai_question"
    )
    
    if st.button("Ask Claude", type="primary", use_container_width=True, key="ask_ai"):
        if not user_question.strip():
            st.warning("Please enter a question.")
        else:
            with st.spinner("Generating SQL query..."):
                try:
                    # Get schema for context
                    schema = """
                    run(run_id, flowcell_id, device_id, exp_start_time, sample_id)
                    run_qc(run_id, n_reads_total, mean_read_length, mapped_reads_total, mapping_pct)
                    barcode_qc(run_id, barcode, total_reads, mapped_reads_total, mapping_pct, mean_read_length, mapq_mean)
                    alignment_groups(run_id, ref_name, read_count, mapq_mean)
                    alignment_groups_by_barcode(run_id, barcode, ref_name, read_count, mapq_mean)
                    sequence_clusters(run_id, seq_hash, cluster_id, seq_len)
                    cluster_counts(run_id, barcode, cluster_id, count)
                    """
                    
                    generated_sql = generate_sql_from_natural_language(user_question, schema)
                    
                    st.success("Generated SQL query:")
                    st.code(generated_sql, language="sql")
                    
                    # Auto-populate the query text area
                    st.session_state["query_text"] = generated_sql
                    
                    # Auto-run the query
                    try:
                        df = sql_df(conn, generated_sql)
                        if df.empty:
                            st.info("Query returned no rows.")
                        else:
                            st.success(f"{len(df):,} rows returned")
                            st.dataframe(df, use_container_width=True)
                            download_button_df("Download as CSV", df, "ai_query_results.csv")
                    except Exception as e:
                        st.error(f"Query error: {e}")
                        st.info("The SQL was generated but failed to execute. You can edit it manually below.")
                        
                except Exception as e:
                    st.error(f"Failed to generate query: {e}")

    st.divider()
    st.subheader("📝 Manual SQL Query (Advanced)")

    # Schema reference
    with st.expander("📋 Schema reference", expanded=False):
        st.markdown("""
| Table | Key columns |
|---|---|
| `run` | `run_id`, `flowcell_id`, `device_id`, `exp_start_time`, `sequencing_kit`, `sample_id` |
| `run_qc` | `run_id`, `n_reads_total`, `mean_read_length`, `n_primary`, `mapped_reads_total`, `mapping_pct` |
| `barcode_qc` | `run_id`, `barcode`, `total_reads`, `mapped_reads_total`, `mapping_pct`, `mean_read_length`, `mapq_mean` |
| `alignment_groups` | `run_id`, `ref_name`, `read_count`, `mapq_mean` |
| `alignment_groups_by_barcode` | `run_id`, `barcode`, `ref_name`, `read_count`, `mapq_mean` |
| `sequence_clusters` | `run_id`, `seq_hash`, `cluster_id`, `seq_len` |
| `cluster_counts` | `run_id`, `barcode`, `cluster_id`, `count` |
        """)

    # Example queries
    with st.expander("💡 Example queries", expanded=False):
        examples = {
            "All runs": "SELECT * FROM run;",
            "QC summary across all runs": "SELECT r.run_id, r.sample_id, r.exp_start_time, q.n_reads_total, q.mapping_pct, q.mean_read_length\nFROM run r\nJOIN run_qc q ON r.run_id = q.run_id\nORDER BY r.exp_start_time DESC;",
            "Top barcodes by mapping % for a run": "SELECT barcode, total_reads, mapped_reads_total, mapping_pct, mean_read_length\nFROM barcode_qc\nWHERE run_id = 'YOUR_RUN_ID'\nORDER BY mapping_pct DESC\nLIMIT 20;",
            "Top clusters in a specific barcode": "SELECT cluster_id, count FROM cluster_counts WHERE run_id = 'YOUR_RUN_ID' AND barcode = 'barcode08' ORDER BY count DESC LIMIT 10;",
            "Barcodes with low mapping %": "SELECT run_id, barcode, total_reads, mapping_pct\nFROM barcode_qc\nWHERE mapping_pct < 50\nORDER BY mapping_pct ASC;",
        }
        for label, sql in examples.items():
            if st.button(label, use_container_width=True, key=f"example_{label}"):
                st.session_state["query_text"] = sql

    # Query input
    query = st.text_area(
        "SQL query",
        value=st.session_state.get("query_text", "SELECT * FROM run;"),
        height=150,
        key="query_text",
    )

    if st.button("Run query", type="primary", use_container_width=True):
        if not query.strip():
            st.warning("Enter a query first.")
            return

        # Only allow read-only queries
        first_word = query.strip().split()[0].upper()
        if first_word not in ("SELECT", "WITH", "EXPLAIN", "PRAGMA"):
            st.error("Only SELECT queries are allowed. This interface is read-only.")
            return

        try:
            df = sql_df(conn, query)
            if df.empty:
                st.info("Query returned no rows.")
            else:
                st.success(f"{len(df):,} rows returned")
                st.dataframe(df, use_container_width=True)
                download_button_df("Download as CSV", df, "query_results.csv")

        except Exception as e:
            st.error(f"Query error: {e}")


def page_candidates(conn: sqlite3.Connection, db_path: Path):
    st.header("⭐ Top Abundant Clusters")
    st.caption("Ranks clusters by read count within each barcode. Use this to identify the most abundant sequences in each round/sample.")

    if sql_df(conn, "SELECT 1 FROM run LIMIT 1;").empty:
        st.info("No runs in the database yet.")
        return

    run_id = run_selectbox(conn, "Select run", key="candidates_run_select")
    if run_id is None:
        return

    # Check if clustering has been run
    any_clusters = sql_df(conn, "SELECT 1 FROM cluster_counts WHERE run_id=? LIMIT 1;", (run_id,))

    # Show clustering form (always, whether or not clustering exists)
    with st.expander("🔧 Clustering Settings" + (" (re-run to change parameters)" if not any_clusters.empty else ""), expanded=any_clusters.empty):
        # Show system RAM info
        total_ram = get_total_ram_gb()
        available_ram = get_available_ram_gb()
        recommended_cap = recommend_max_sequences()
        
        st.info(f"""
        **System RAM:** {total_ram:.1f} GB total, {available_ram:.1f} GB available  
        **Recommended clustering cap:** {recommended_cap:,} sequences  
        _(More RAM = more sequences = better chance of finding enrichment)_
        """)

        zip_path_disk = st.text_input(
            "Zip path on disk",
            value=st.session_state.get("last_zip_path", ""),
            placeholder="/Users/you/Downloads/run.zip",
            key="cluster_zip_path"
        )
        work_dir_str = st.text_input(
            "Work dir",
            value=st.session_state.get("last_work_dir", str(Path(tempfile.gettempdir()) / "minion_tmp")),
            key="cluster_work_dir"
        )
        
        col1, col2, col3 = st.columns(3)
        with col1:
            min_len = st.number_input("Min sequence length (bp)", min_value=100, value=MIN_SEQ_LEN, step=10, key="cluster_min_len")
            min_mapq = st.number_input("Min MAPQ", min_value=0, max_value=60, value=10, step=1, key="cluster_min_mapq")
        with col2:
            max_len = st.number_input("Max sequence length (bp)", min_value=100, value=MAX_SEQ_LEN, step=10, key="cluster_max_len")
            min_seq_id = st.number_input("MMseqs2 min sequence identity", min_value=0.5, max_value=1.0, value=0.80, step=0.01, 
                                          help="Proteins cluster at 80% by default. Lower = more permissive clustering.", key="cluster_min_seq_id")
        with col3:
            max_seqs_cap = st.number_input("Max sequences to cluster", min_value=1000, max_value=200_000, value=recommended_cap, step=5000,
                                            help=f"Recommended: {recommended_cap:,} based on your RAM. Higher values may cause out-of-memory errors.", key="cluster_max_seqs")

        button_label = "Re-run clustering" if not any_clusters.empty else "Run clustering"
        if st.button(button_label, type="primary", use_container_width=True, key="cluster_run_button"):
            zp = Path(zip_path_disk.strip())
            if not zp.exists():
                st.error("Zip path does not exist.")
                return

            check_mmseqs()

            work_dir = Path(work_dir_str)
            work_dir.mkdir(parents=True, exist_ok=True)

            entries = find_bams_in_zip(zp)
            if not entries:
                st.error("No BAMs found in zip.")
                return

            # Delete old clustering if re-running
            if not any_clusters.empty:
                st.warning("Deleting old clustering results...")
                conn.execute("DELETE FROM cluster_counts WHERE run_id = ?", (run_id,))
                conn.execute("DELETE FROM sequence_clusters WHERE run_id = ?", (run_id,))
                conn.commit()

            status = st.empty()
            progress = st.progress(0)

            bam_paths = []
            with zipfile.ZipFile(zp, "r") as z:
                for i, entry in enumerate(entries, start=1):
                    tmp_dir = work_dir / f"tmp_cluster_{i:03d}"
                    if tmp_dir.exists():
                        shutil.rmtree(tmp_dir, ignore_errors=True)
                    tmp_dir.mkdir(parents=True, exist_ok=True)
                    bam_path = safe_extract_member(z, entry.bam_member, tmp_dir)
                    if entry.bai_member:
                        safe_extract_member(z, entry.bai_member, tmp_dir)
                    bam_paths.append((bam_path, entry.barcode))

            # Pass custom parameters to clustering
            build_sequence_clusters(
                run_id=run_id,
                bam_paths=bam_paths,
                conn=conn,
                work_dir=work_dir,
                min_mapq=int(min_mapq),
                min_len=int(min_len),
                max_len=int(max_len),
                max_seqs=int(max_seqs_cap),
                ui_status=status,
                ui_progress=progress,
            )

            for bam_path, _ in bam_paths:
                shutil.rmtree(bam_path.parent, ignore_errors=True)

            st.rerun()

    # If no clustering exists, stop here
    if any_clusters.empty:
        return

    # Clustering exists - show abundance ranking
    barcodes_df = sql_df(conn, """
        SELECT barcode, total_reads
        FROM barcode_qc
        WHERE run_id = ?
        ORDER BY barcode;
    """, (run_id,))

    counts_df = sql_df(conn, """
        SELECT barcode, cluster_id, count
        FROM cluster_counts
        WHERE run_id = ?
        ORDER BY barcode, count DESC;
    """, (run_id,))

    # Build abundance table
    top_per_barcode = []
    for _, row in barcodes_df.iterrows():
        bc = row["barcode"]
        top_clusters = counts_df[counts_df["barcode"] == bc].nlargest(100, "count")  # Top 100 per barcode
        for rank, (_, cluster_row) in enumerate(top_clusters.iterrows(), start=1):
            top_per_barcode.append({
                "barcode": bc,
                "rank": rank,
                "cluster_id": cluster_row["cluster_id"],
                "reads": int(cluster_row["count"]),
                "pct_of_barcode": round(cluster_row["count"] / row["total_reads"] * 100, 3) if row["total_reads"] > 0 else 0,
            })

    if not top_per_barcode:
        st.error("No abundance data to display.")
        return

    abundance_df = pd.DataFrame(top_per_barcode)

    # Filters
    st.subheader("Filters")
    col1, col2 = st.columns(2)
    with col1:
        min_reads = st.number_input("Min reads", min_value=0, value=10, step=5)
    with col2:
        max_rank = st.number_input("Max rank (per barcode)", min_value=1, max_value=100, value=20, step=5,
                                    help="Show top N clusters per barcode")

    filtered = abundance_df[
        (abundance_df["reads"] >= min_reads) &
        (abundance_df["rank"] <= max_rank)
    ]

    st.subheader(f"{len(filtered)} clusters (top {max_rank} per barcode, ≥{min_reads} reads)")
    st.dataframe(filtered, use_container_width=True)
    download_button_df("Download abundant clusters CSV", filtered, f"{run_id}_abundant_clusters.csv")

    # Barcode summary stats
    with st.expander("📊 Barcode summary statistics"):
        summary = []
        for bc in barcodes_df["barcode"]:
            bc_clusters = counts_df[counts_df["barcode"] == bc]
            summary.append({
                "barcode": bc,
                "n_clusters": len(bc_clusters),
                "total_reads_clustered": int(bc_clusters["count"].sum()),
                "max_cluster_size": int(bc_clusters["count"].max()) if not bc_clusters.empty else 0,
                "mean_cluster_size": round(bc_clusters["count"].mean(), 1) if not bc_clusters.empty else 0,
            })
        summary_df = pd.DataFrame(summary)
        st.dataframe(summary_df, use_container_width=True)
def main():
    st.set_page_config(page_title=APP_TITLE, layout="wide")
    st.title(APP_TITLE)

    st.sidebar.header("Database")
    db_path_str = st.sidebar.text_input("SQLite DB path", value=DEFAULT_DB)
    db_path = Path(db_path_str).expanduser().resolve()

    conn = connect_db(db_path)

    st.sidebar.divider()
    page = st.sidebar.radio(
        "Navigate",
        ["Overview", "Run detail", "Candidates", "Query", "Ingest"],
        index=0,
    )

    if page == "Overview":
        page_overview(conn)
    elif page == "Run detail":
        page_run_detail(conn)
    elif page == "Candidates":
        page_candidates(conn, db_path)
    elif page == "Query":
        page_query(conn)
    elif page == "Ingest":
        page_ingest(conn, db_path)

    st.sidebar.divider()
    st.sidebar.caption("Tip: Keep one 'master' DB (e.g., minion_master.db) and ingest all runs into it.")


if __name__ == "__main__":
    main()
