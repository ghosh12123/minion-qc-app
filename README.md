# MinION Nanobody Analysis v2

A Streamlit web app for analyzing MinION nanobody sequencing data. Ingests FASTQ run folders, runs the trim-translate-cluster-normalize pipeline, and provides interactive enrichment analysis across panning rounds.

## Requirements

- Python 3.10 or higher
- MMseqs2 (see below)

## Installation

### 1. Install Python dependencies

```bash
pip install -r requirements.txt
```

### 2. Install MMseqs2

MMseqs2 is the sequence clustering tool used by the pipeline. It must be installed separately.

**Mac (Homebrew):**
```bash
brew install mmseqs2
```

If you don't have Homebrew installed:
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

**Linux:**
```bash
conda install -c conda-forge -c bioconda mmseqs2
```
Or download a pre-built binary from the [MMseqs2 releases page](https://github.com/soedinglab/MMseqs2/releases).

**Windows (WSL):**
```bash
conda install -c conda-forge -c bioconda mmseqs2
```

**Verify installation:**
```bash
mmseqs version
```

## Running the App

```bash
streamlit run fastq_analysis_v2.py
```

Then open your browser and go to `http://localhost:8501`.

## Usage

1. Go to the **Ingest** page and enter the path to your run folder (e.g. `/Users/you/Downloads/run5`). The app auto-detects all barcode folders and groups them by target.
2. Set your START and END anchor sequences, minimum AA length, and MMseqs2 parameters.
3. If your run has no TG1 barcode of its own, use the **External TG1** section to select a TG1 from an already-ingested run or provide a folder path directly.
4. Click **Start Ingest**. Results are stored in a local SQLite database (`~/minion_nanobody_v2.db`) and never need to be recomputed.
5. Go to the **Enrichment** page, select a run and target, and explore the results.

## Notes

- Default anchor sequences (`ATGGCC` / `GGCGCGC`) are construct-specific. Change them if your library uses a different vector.
- `easy-cluster` (default) matches the lab's existing pipeline exactly. `easy-linclust` is ~10x faster but slightly less accurate.
- The baseline threshold should be adjusted based on TG1 sequencing depth. Lower TG1 depth = lower threshold.
