# Bioinformatics Sequence Alignment

A hybrid DNA similarity search project using NCBI BLASTN remotely or BLAST+ locally.

## Fast default mode

Run:

```bash
python main.py
```

Paste a DNA sequence and press Enter twice. The program submits the query to NCBI and immediately prints a request ID and result link. It does not block while NCBI processes the search.

After the result is ready, retrieve parsed hits with:

```bash
python main.py --check RID_HERE
```

The program will ask for the original query length to calculate coverage.

## Local fast mode

With BLAST+ installed and a local database:

```bash
python main.py --db data/nt
```

Local mode is recommended for repeated searches. NCBI Datasets can download targeted gene or genome data, which can then be formatted as a local BLAST database.

The remote mode uses NCBI's asynchronous BLAST URL API. The local mode avoids the public queue entirely.
