"""Interactive and command-line DNA similarity search."""
import argparse
from pathlib import Path
import sys

from src.fasta import read_fasta
from src.ncbi import NCBIClient
from src.search import search_references


def normalize_sequence(raw: str) -> str:
    sequence = "".join(raw.split()).upper()
    if not sequence:
        raise ValueError("Sequence cannot be empty.")
    invalid = set(sequence) - set("ACGTN")
    if invalid:
        raise ValueError(f"Invalid DNA symbols: {sorted(invalid)}")
    return sequence


def interactive() -> int:
    print("DNA Sequence Alignment - interactive mode")
    email = input("NCBI contact email: ").strip()
    gene = input("Gene name to check in NCBI (for example BRCA1): ").strip()
    sequence = normalize_sequence(input("Paste the DNA sequence: "))
    organism = input("Organism filter (optional, for example Homo sapiens): ").strip()
    client = NCBIClient(email=email)
    result = client.gene_exists(gene, organism=organism)
    if not result.exists:
        print(f"NCBI Gene check: no matching gene found for '{gene}'.")
        return 0
    print(f"NCBI Gene check: found {result.count} record(s) for '{gene}'.")
    print("Gene IDs:", ", ".join(result.ids[:10]))
    term = f"{gene}[Gene Name]"
    if organism:
        term += f' AND "{organism}"[Organism]'
    print("Searching NCBI nucleotide records related to the gene...")
    references = client.search_and_fetch(term, retmax=10)
    print(f"Downloaded {len(references)} reference sequence(s).")
    hits = search_references(sequence, references, k=min(11, len(sequence)), min_seeds=1)
    for rank, hit in enumerate(hits[:10], 1):
        print(f"{rank}. {hit.name} | score={hit.score} | identity={hit.identity:.2f}% | seeds={hit.seeds}")
        print(f"   {hit.aligned_query}")
        print(f"   {hit.aligned_reference}")
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description="Search local FASTA files or NCBI nucleotide records.")
    p.add_argument("--interactive", action="store_true", help="Prompt for a sequence and NCBI gene search")
    source = p.add_mutually_exclusive_group()
    source.add_argument("--reference", type=Path)
    source.add_argument("--ncbi", metavar="TERM")
    p.add_argument("--query", type=Path)
    p.add_argument("--sequence", help="DNA sequence supplied directly on the command line")
    p.add_argument("--limit", type=int, default=10)
    p.add_argument("--kmer", type=int, default=11)
    p.add_argument("--min-seeds", type=int, default=1)
    p.add_argument("--email")
    p.add_argument("--api-key")
    args = p.parse_args()
    if args.interactive:
        return interactive()
    if not args.query and not args.sequence:
        p.error("provide --interactive, --query, or --sequence")
    if not args.reference and not args.ncbi:
        p.error("provide --reference or --ncbi")
    try:
        queries = [("typed_sequence", normalize_sequence(args.sequence))] if args.sequence else read_fasta(args.query)
        if args.reference:
            references = read_fasta(args.reference)
        else:
            if not args.email:
                raise ValueError("--email is required with --ncbi")
            references = NCBIClient(args.email, args.api_key).search_and_fetch(args.ncbi, args.limit)
        for query_name, query in queries:
            hits = search_references(query, references, k=min(args.kmer, len(query)), min_seeds=args.min_seeds)
            print(f"\nQuery: {query_name} ({len(query)} bp)")
            for rank, hit in enumerate(hits, 1):
                print(f"{rank}. {hit.name} | score={hit.score} | identity={hit.identity:.2f}% | seeds={hit.seeds}")
        return 0
    except (OSError, ValueError, RuntimeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
