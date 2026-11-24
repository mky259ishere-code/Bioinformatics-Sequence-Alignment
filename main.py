"""Fast local DNA similarity search using NCBI BLAST+."""
import argparse
import sys

from src.fasta import read_fasta
from src.local_blast import run_local_blast


def clean_sequence(raw: str) -> str:
    sequence = "".join(raw.split()).upper()
    if not sequence:
        raise ValueError("The DNA sequence cannot be empty.")
    invalid = set(sequence) - set("ACGTN")
    if invalid:
        raise ValueError(f"Invalid DNA symbols: {sorted(invalid)}")
    return sequence


def get_sequence(args) -> str:
    if args.sequence:
        return clean_sequence(args.sequence)
    if args.fasta:
        records = read_fasta(args.fasta)
        if len(records) != 1:
            raise ValueError("Please provide exactly one query sequence.")
        return records[0][1]
    print("Paste a DNA sequence. Press Enter twice when finished:")
    lines = []
    while True:
        line = input()
        if not line.strip():
            break
        lines.append(line)
    return clean_sequence("".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description="Find similar DNA sequences using local BLASTN.")
    parser.add_argument("--db", required=True, help="Local BLAST database prefix, e.g. data/nt")
    parser.add_argument("--sequence")
    parser.add_argument("--fasta")
    parser.add_argument("--max-hits", type=int, default=20)
    parser.add_argument("--blastn", default="blastn", help="blastn executable or full path")
    args = parser.parse_args()
    try:
        sequence = get_sequence(args)
        print(f"Searching local BLAST database with {len(sequence)} bases...")
        hits = run_local_blast(sequence, args.db, args.max_hits, args.blastn)
        if not hits:
            print("No significant matches were returned.")
            return 0
        print(f"\nFound {len(hits)} possible matches:\n")
        for index, hit in enumerate(hits, 1):
            print(f"{index}. {hit.title}")
            print(f"   Accession: {hit.accession}")
            print(f"   Identity: {hit.identity:.2f}%")
            print(f"   Query coverage: {hit.coverage:.2f}%")
            print(f"   Alignment length: {hit.alignment_length}")
            print(f"   E-value: {hit.evalue}")
            print(f"   NCBI link: {hit.link}\n")
        return 0
    except (OSError, ValueError, RuntimeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
