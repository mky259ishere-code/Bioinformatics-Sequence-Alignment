"""NCBI-backed, BLAST-like DNA similarity search."""
import argparse
from pathlib import Path
import sys

from src.fasta import read_fasta
from src.ncbi import NCBIClient
from src.search import search_references


def main() -> int:
    p = argparse.ArgumentParser(description="Search local FASTA files or NCBI nucleotide records.")
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--reference", type=Path, help="Local reference FASTA file")
    source.add_argument("--ncbi", metavar="TERM", help="NCBI nucleotide search term, e.g. 'BRCA1[Gene] AND Homo sapiens[Organism]'")
    p.add_argument("--query", type=Path, required=True, help="Query FASTA file")
    p.add_argument("--limit", type=int, default=10, help="Maximum NCBI records to retrieve")
    p.add_argument("--kmer", type=int, default=11, help="Seed length")
    p.add_argument("--min-seeds", type=int, default=1)
    p.add_argument("--email", help="Contact email required/recommended by NCBI for automated use")
    p.add_argument("--api-key", help="Optional NCBI API key")
    args = p.parse_args()

    try:
        queries = read_fasta(args.query)
        if args.reference:
            references = read_fasta(args.reference)
        else:
            if not args.email:
                raise ValueError("--email is required when using --ncbi")
            client = NCBIClient(email=args.email, api_key=args.api_key)
            references = client.search_and_fetch(args.ncbi, retmax=args.limit)
        for query_name, query in queries:
            hits = search_references(query, references, k=args.kmer, min_seeds=args.min_seeds)
            print(f"\nQuery: {query_name} ({len(query)} bp)")
            if not hits:
                print("No seed-supported hits found.")
                continue
            print("Rank  Reference                         Score  Identity  Seeds  Span")
            for rank, hit in enumerate(hits, 1):
                print(f"{rank:<5} {hit.name[:33]:<33} {hit.score:>6}  {hit.identity:>7.2f}%  {hit.seeds:>5}  {hit.span:>4}")
                print(f"      {hit.aligned_query}")
                print(f"      {hit.aligned_reference}")
        return 0
    except (OSError, ValueError, RuntimeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
