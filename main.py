"""User-friendly NCBI BLASTN sequence search."""
import argparse
import sys

from src.blast import BlastClient
from src.fasta import read_fasta


def clean_sequence(raw: str) -> str:
    sequence = "".join(raw.split()).upper()
    if not sequence:
        raise ValueError("The DNA sequence cannot be empty.")
    invalid = set(sequence) - set("ACGTN")
    if invalid:
        raise ValueError(f"Invalid DNA symbols: {sorted(invalid)}")
    return sequence


def main() -> int:
    parser = argparse.ArgumentParser(description="Find similar DNA sequences in NCBI using BLASTN.")
    parser.add_argument("--sequence", help="DNA sequence to search")
    parser.add_argument("--fasta", help="FASTA file containing the query sequence")
    parser.add_argument("--max-hits", type=int, default=20)
    args = parser.parse_args()

    try:
        if args.sequence:
            sequence = clean_sequence(args.sequence)
        elif args.fasta:
            records = read_fasta(args.fasta)
            if len(records) != 1:
                raise ValueError("Please provide exactly one query sequence.")
            sequence = records[0][1]
        else:
            print("Paste a DNA sequence. Press Enter twice when finished:")
            lines = []
            while True:
                line = input()
                if not line.strip():
                    break
                lines.append(line)
            sequence = clean_sequence("".join(lines))

        print(f"Submitting {len(sequence)} bases to NCBI BLASTN...")
        hits = BlastClient().search(sequence, max_hits=args.max_hits)
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
    except (OSError, ValueError, RuntimeError, TimeoutError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
