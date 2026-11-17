"""Command-line DNA similarity search using local alignment."""

import argparse
from pathlib import Path
import sys

from src.alignment import smith_waterman
from src.fasta import read_fasta


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare a query DNA sequence against reference FASTA sequences.")
    parser.add_argument("query", type=Path, help="FASTA file containing one or more query sequences")
    parser.add_argument("reference", type=Path, help="FASTA file containing reference sequences")
    parser.add_argument("--match", type=int, default=2)
    parser.add_argument("--mismatch", type=int, default=-1)
    parser.add_argument("--gap", type=int, default=-2)
    args = parser.parse_args()

    try:
        queries = read_fasta(args.query)
        references = read_fasta(args.reference)
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    for query_name, query_sequence in queries:
        results = [(reference_name, smith_waterman(query_sequence, reference_sequence, args.match, args.mismatch, args.gap)) for reference_name, reference_sequence in references]
        results.sort(key=lambda item: (item[1].score, item[1].identity_percent), reverse=True)
        print(f"\\nQuery: {query_name} ({len(query_sequence)} bp)")
        print("Rank  Reference                 Score  Identity  Alignment")
        print("----  ------------------------  -----  --------  ---------")
        for rank, (reference_name, result) in enumerate(results, 1):
            print(f"{rank:<5} {reference_name[:24]:<24}  {result.score:>5}  {result.identity_percent:>7.2f}%  {result.alignment_length:>9}")
            print(f"      {result.aligned_query}")
            print(f"      {result.aligned_reference}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
