"""Streaming FASTA parsing with validation."""
from pathlib import Path
from typing import Iterator


def read_fasta_iter(path: str | Path) -> Iterator[tuple[str, str]]:
    header = None
    parts: list[str] = []
    for line_number, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(parts)
            header, parts = line[1:].strip() or f"sequence_{line_number}", []
        elif header is None:
            raise ValueError(f"Line {line_number}: sequence before FASTA header")
        else:
            parts.append("".join(line.split()).upper())
    if header is not None:
        yield header, "".join(parts)


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    records = list(read_fasta_iter(path))
    if not records:
        raise ValueError("FASTA file contains no sequences")
    for name, sequence in records:
        if not sequence:
            raise ValueError(f"Sequence '{name}' is empty")
        invalid = set(sequence) - set("ACGTN")
        if invalid:
            raise ValueError(f"Sequence '{name}' contains invalid bases: {sorted(invalid)}")
    return records
