#!/usr/bin/env python3
"""Validate the headerless three-column prediction contract used by LAFA."""

import argparse
import gzip
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lafa_main import read_fasta


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("prediction_file")
    parser.add_argument("--query_file", help="Optional FASTA used to verify query IDs")
    args = parser.parse_args()
    path = Path(args.prediction_file)
    opener = gzip.open if path.suffix == ".gz" else open
    rows = 0
    proteins = set()
    pairs = set()
    with opener(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 3:
                raise SystemExit(f"line {line_number}: expected 3 tab-separated columns")
            protein, term, score_text = fields
            if not protein or not term.startswith("GO:") or len(term) != 10 or not term[3:].isdigit():
                raise SystemExit(f"line {line_number}: invalid protein or GO term")
            try:
                score = float(score_text)
            except ValueError:
                raise SystemExit(f"line {line_number}: invalid score {score_text!r}") from None
            if not math.isfinite(score) or not 0 <= score <= 1:
                raise SystemExit(f"line {line_number}: score outside [0, 1]")
            if (protein, term) in pairs:
                raise SystemExit(f"line {line_number}: duplicate pair {protein}, {term}")
            pairs.add((protein, term))
            proteins.add(protein)
            rows += 1
    coverage = ""
    if args.query_file:
        query_ids = {protein for protein, _ in read_fasta(Path(args.query_file))}
        unknown = proteins - query_ids
        if unknown:
            examples = ", ".join(sorted(unknown)[:5])
            raise SystemExit(f"prediction IDs absent from query FASTA: {examples}")
        coverage = f", {len(query_ids - proteins)} query proteins without predictions"
    print(f"OK: {rows} rows, {len(proteins)} proteins, headerless 3-column TSV{coverage}")


if __name__ == "__main__":
    main()
