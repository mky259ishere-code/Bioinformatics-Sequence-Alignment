# Bioinformatics Sequence Alignment

A teaching project for DNA similarity search.

## Current workflow

Run the program without arguments:

```bash
python main.py
```

Paste a DNA sequence and press Enter twice. The program submits it to NCBI BLASTN and reports ranked possible matches with accession links, identity, query coverage, alignment length, and E-value.

You can also use a command-line sequence:

```bash
python main.py --sequence ATGCGTACGTTAGCTAGC
```

Or a one-record FASTA file:

```bash
python main.py --fasta data/sample_query.fasta
```

The program uses NCBI's public BLAST URL API. Searches may take time because NCBI processes them remotely. It is intended for individual educational use, not high-volume batch searches. NCBI recommends respecting request limits and identifying automated tools; this version does not ask the user for an email.

This project is an educational client and does not replace the official NCBI BLAST website.
