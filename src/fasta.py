"""Utilities for reading FASTA files."""

from pathlib import Path


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Read a FASTA file and return (header, sequence) records."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence_parts: list[str] = []

    for line_number, raw_line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence_parts)))
            header = line[1:].strip() or f"sequence_{len(records) + 1}"
            sequence_parts = []
        elif header is None:
            raise ValueError(f"Line {line_number}: sequence data appears before a FASTA header")
        else:
            sequence_parts.append(line.replace(" ", "").upper())

    if header is not None:
        records.append((header, "".join(sequence_parts)))

    if not records:
        raise ValueError("The FASTA file contains no sequences")

    for name, sequence in records:
        if not sequence:
            raise ValueError(f"Sequence '{name}' is empty")
        invalid = set(sequence) - set("ACGTN")
        if invalid:
            raise ValueError(f"Sequence '{name}' contains invalid bases: {sorted(invalid)}")

    return records
