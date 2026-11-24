"""Run BLAST+ locally and parse tabular output."""
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote


@dataclass(frozen=True)
class LocalBlastHit:
    accession: str
    title: str
    identity: float
    coverage: float
    alignment_length: int
    evalue: str
    bitscore: float

    @property
    def link(self) -> str:
        return f"https://www.ncbi.nlm.nih.gov/nuccore/{quote(self.accession)}"


def run_local_blast(sequence: str, database: str, max_hits: int = 20,
                    executable: str = "blastn") -> list[LocalBlastHit]:
    command = shutil.which(executable) or executable
    fields = "qseqid saccver stitle pident length qlen qstart qend evalue bitscore"
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False, encoding="utf-8") as query:
        query.write(">user_query\n" + sequence + "\n")
        query_path = query.name
    try:
        result = subprocess.run(
            [command, "-query", query_path, "-db", database, "-task", "megablast",
             "-max_target_seqs", str(max_hits), "-outfmt", f"6 {fields}"],
            capture_output=True, text=True, check=False,
        )
    except FileNotFoundError as error:
        raise RuntimeError(
            "blastn was not found. Install NCBI BLAST+ and add its bin folder to PATH."
        ) from error
    finally:
        Path(query_path).unlink(missing_ok=True)
    if result.returncode:
        message = result.stderr.strip() or "blastn failed."
        raise RuntimeError(message)
    hits = []
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != 10:
            continue
        _, accession, title, identity, length, qlen, qstart, qend, evalue, bitscore = parts
        query_span = abs(int(qend) - int(qstart)) + 1
        hits.append(LocalBlastHit(
            accession, title, float(identity), 100.0 * query_span / int(qlen),
            int(length), evalue, float(bitscore),
        ))
    return hits
