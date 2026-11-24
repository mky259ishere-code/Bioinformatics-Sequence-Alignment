"""Hybrid DNA similarity search: asynchronous NCBI BLAST or local BLAST+."""
import argparse
import sys

from src.blast import BlastClient
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


def print_hits(hits) -> None:
    if not hits:
        print("No significant matches were returned.")
        return
    print(f"\nFound {len(hits)} possible matches:\n")
    for index, hit in enumerate(hits, 1):
        print(f"{index}. {hit.title}")
        print(f"   Accession: {hit.accession}")
        print(f"   Identity: {hit.identity:.2f}%")
        print(f"   Query coverage: {hit.coverage:.2f}%")
        print(f"   Alignment length: {hit.alignment_length}")
        print(f"   E-value: {hit.evalue}")
        print(f"   NCBI link: {hit.link}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Find similar DNA sequences in NCBI or a local BLAST database.")
    parser.add_argument("--db", help="Optional local BLAST database prefix")
    parser.add_argument("--check", metavar="RID", help="Retrieve completed results for an NCBI request ID")
    parser.add_argument("--sequence")
    parser.add_argument("--fasta")
    parser.add_argument("--max-hits", type=int, default=20)
    parser.add_argument("--blastn", default="blastn")
    args = parser.parse_args()
    try:
        if args.check:
            length = int(input("Enter the original query length: "))
            print_hits(BlastClient().parse_results(BlastClient().fetch_xml(args.check), length))
            return 0
        sequence = get_sequence(args)
        if args.db:
            print(f"Searching local BLAST database with {len(sequence)} bases...")
            print_hits(run_local_blast(sequence, args.db, args.max_hits, args.blastn))
        else:
            rid = BlastClient().submit(sequence, args.max_hits)
            print("\nNCBI BLAST search submitted successfully.")
            print(f"Request ID: {rid}")
            print(f"Open this link when ready: {BlastClient.RESULT_URL.format(rid)}")
            print(f"To retrieve parsed results later: python main.py --check {rid}")
        return 0
    except (OSError, ValueError, RuntimeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
