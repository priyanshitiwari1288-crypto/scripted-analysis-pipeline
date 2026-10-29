#!/usr/bin/env python3
"""
Small, reproducible DNA FASTA analysis pipeline.
Run from the project root: python pipeline.py
Only uses the Python standard library.
"""

# --------------------------- PARAMETERS ---------------------------
INPUT_FASTA = "data/sample.fasta"
OUTPUT_DIR = "results"
KMER_SIZE = 3
# -----------------------------------------------------------------

from pathlib import Path
import csv
import json
import re
import sys

DNA_ALPHABET = set("ACGTN")

def read_fasta(path):
    """Read a FASTA file and return records as (header, sequence) tuples."""
    records = []
    header = None
    chunks = []
    with open(path, "r", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks).upper()))
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header on line {line_no}")
                chunks = []
            else:
                if header is None:
                    raise ValueError(
                        f"Sequence data appears before a FASTA header (line {line_no})"
                    )
                chunks.append(re.sub(r"\s+", "", line))
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    if not records:
        raise ValueError("No FASTA records found.")
    for name, seq in records:
        if not seq:
            raise ValueError(f"Sequence '{name}' is empty.")
        invalid = sorted(set(seq) - DNA_ALPHABET)
        if invalid:
            raise ValueError(
                f"Sequence '{name}' contains unsupported symbols: {''.join(invalid)}. "
                "This starter pipeline expects DNA symbols A, C, G, T, or N."
            )
    return records

def count_bases(sequence):
    return {base: sequence.count(base) for base in "ACGTN"}

def gc_percent(sequence):
    # N bases are excluded from the denominator because they are ambiguous.
    known = sum(sequence.count(base) for base in "ACGT")
    if known == 0:
        return 0.0
    return 100.0 * (sequence.count("G") + sequence.count("C")) / known

def kmer_counts(sequence, k):
    counts = {}
    for i in range(len(sequence) - k + 1):
        word = sequence[i:i+k]
        if "N" in word:
            continue
        counts[word] = counts.get(word, 0) + 1
    return dict(sorted(counts.items()))

def main():
    input_path = Path(INPUT_FASTA)
    out_dir = Path(OUTPUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)

    if KMER_SIZE < 1:
        raise ValueError("KMER_SIZE must be at least 1.")
    if not input_path.exists():
        raise FileNotFoundError(
            f"Input file not found: {input_path}. Add a FASTA file at this path."
        )

    # Step 1: parse and validate the input FASTA.
    records = read_fasta(input_path)

    # Step 2: calculate per-sequence summary statistics.
    summary_rows = []
    all_kmers = {}
    for name, seq in records:
        row = {
            "sequence_id": name,
            "length": len(seq),
            "gc_percent": round(gc_percent(seq), 3),
            **{f"count_{base}": n for base, n in count_bases(seq).items()},
        }
        summary_rows.append(row)
        for word, count in kmer_counts(seq, KMER_SIZE).items():
            all_kmers[word] = all_kmers.get(word, 0) + count

    # Step 3: keep intermediate results for transparency/reproducibility.
    normalized_path = out_dir / "normalized_sequences.fasta"
    with normalized_path.open("w", encoding="utf-8") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                handle.write(seq[i:i+60] + "\n")

    summary_path = out_dir / "sequence_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    kmer_path = out_dir / "kmer_counts.csv"
    with kmer_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["kmer", "count"])
        for word, count in sorted(all_kmers.items(), key=lambda item: (-item[1], item[0])):
            writer.writerow([word, count])

    run_info = {
        "input_file": str(input_path),
        "sequence_count": len(records),
        "kmer_size": KMER_SIZE,
        "outputs": [
            str(normalized_path), str(summary_path), str(kmer_path)
        ],
    }
    with (out_dir / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(run_info, handle, indent=2)

    print("Pipeline completed successfully.")
    print(f"Sequences analysed: {len(records)}")
    print(f"Results saved in: {out_dir.resolve()}")
    for output in run_info["outputs"]:
        print(f" - {output}")

if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
