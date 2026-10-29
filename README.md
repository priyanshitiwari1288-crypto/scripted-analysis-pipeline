# Task 03 — Scripted Analysis Pipeline

## Overview
This project turns a small multi-step DNA sequence analysis into a reproducible Python pipeline that runs from one command.

## Requirements
- Python 3.9 or newer
- No third-party packages required (Python standard library only)

## Run
From the project root directory:

```bash
python pipeline.py
```

On some systems, use `python3 pipeline.py`.

## Input
The default input is `data/sample.fasta`. It contains three **illustrative example sequences** so the pipeline can be tested immediately. Replace this file with your own DNA FASTA records when required. The parser expects DNA symbols A, C, G, T, or N.

## Pipeline steps
1. Read and validate the FASTA records.
2. Calculate sequence length, GC percentage, and base counts.
3. Write normalized sequences to an intermediate FASTA file.
4. Write per-sequence statistics to a CSV file.
5. Count k-mers of the configured size and write the results to CSV.
6. Record run settings and output paths in a JSON metadata file.

## Parameters
Edit the constants at the top of `pipeline.py`:
- `INPUT_FASTA`: input FASTA path
- `OUTPUT_DIR`: output directory
- `KMER_SIZE`: k-mer length (default: 3)

## Outputs
Generated in `results/`:
- `normalized_sequences.fasta` — intermediate normalized FASTA
- `sequence_summary.csv` — length, GC percentage, and A/C/G/T/N counts
- `kmer_counts.csv` — k-mer counts across all input sequences
- `run_metadata.json` — run settings and output file list

Intermediate and final outputs are kept so each stage can be inspected.

## Reproducibility and design decisions
- All configurable values are grouped at the top of the script.
- The pipeline uses deterministic processing and Python's standard library.
- Input validation stops the run with a clear error if the FASTA is missing or contains unsupported symbols.
- `N` is treated as an ambiguous base and excluded from the GC-percentage denominator; k-mers containing `N` are skipped.
- The included FASTA data is synthetic/example data for testing, not a biological reference dataset.


