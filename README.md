# Bioinformatics Sequence Alignment

A teaching project for DNA similarity search.

## Features

- Reads ordinary multi-record FASTA files.
- Uses a BLAST-like seed-and-extend workflow: k-mer seeding, seed counting, and Smith-Waterman local extension.
- Searches NCBI's Nucleotide database through E-utilities (ESearch + EFetch).
- Ranks hits by local alignment score, identity, and seed support.
- Uses only Python's standard library.

## Usage

Offline/local search:

```bash
python main.py --query data/sample_query.fasta --reference data/sample_reference.fasta
```

NCBI search:

```bash
python main.py --query data/sample_query.fasta --ncbi "BRCA1[Gene] AND Homo sapiens[Organism]" --limit 10 --email your.email@example.com
```

Use an NCBI API key with `--api-key` if you need more than the default request rate. Do not commit personal API keys to GitHub.

This is an educational BLAST-like implementation, not a replacement for NCBI BLAST. Large-scale searches should use NCBI BLAST or the NCBI Datasets tools.
