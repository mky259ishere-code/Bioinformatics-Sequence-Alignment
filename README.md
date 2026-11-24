# Bioinformatics Sequence Alignment

A teaching project for fast DNA similarity search using local NCBI BLAST+.

## Usage

Install NCBI BLAST+ for Windows from the official NCBI BLAST downloads page, then ensure `blastn` is on PATH.

Create or download a BLAST nucleotide database. A database consists of several files sharing one prefix, such as:

```text
data/nt.nhr
data/nt.nin
data/nt.nsq
```

Run an interactive search:

```bash
python main.py --db data/nt
```

Or provide a sequence directly:

```bash
python main.py --db data/nt --sequence ATGCGTACGTTAGCTAGC
```

The program runs local megablast, then reports ranked matches with accession, identity, query coverage, E-value, and clickable NCBI links. It no longer waits for the public NCBI BLAST queue.

For a small test database, create a FASTA file and format it with:

```bash
makeblastdb -in data/sample_reference.fasta -dbtype nucl -out data/sample_db -parse_seqids
```

NCBI BLAST+ downloads: https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
