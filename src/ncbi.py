"""Rate-limited NCBI Entrez client using ESearch and EFetch."""
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass


@dataclass(frozen=True)
class GeneSearchResult:
    exists: bool
    count: int
    ids: list[str]


class NCBIClient:
    BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    def __init__(self, email: str, api_key: str | None = None, tool: str = "BioinformaticsSequenceAlignment"):
        if not email or "@" not in email:
            raise ValueError("A valid email is required for NCBI requests.")
        self.email, self.api_key, self.tool = email, api_key, tool
        self.last_request = 0.0

    def _get(self, endpoint: str, **params) -> str:
        wait = 0.34 - (time.monotonic() - self.last_request)
        if wait > 0:
            time.sleep(wait)
        params.update(email=self.email, tool=self.tool)
        if self.api_key:
            params["api_key"] = self.api_key
        url = self.BASE + endpoint + "?" + urllib.parse.urlencode(params)
        request = urllib.request.Request(url, headers={"User-Agent": f"{self.tool}/0.1 ({self.email})"})
        with urllib.request.urlopen(request, timeout=30) as response:
            self.last_request = time.monotonic()
            return response.read().decode("utf-8")

    def search_ids(self, database: str, term: str, retmax: int = 10) -> tuple[int, list[str]]:
        xml = ET.fromstring(self._get("esearch.fcgi", db=database, term=term, retmode="xml", retmax=retmax))
        count = int(xml.findtext("./Count", "0"))
        ids = [node.text for node in xml.findall(".//Id") if node.text]
        return count, ids

    def gene_exists(self, gene: str, organism: str = "") -> GeneSearchResult:
        term = f'"{gene}"[Gene Name]'
        if organism:
            term += f' AND "{organism}"[Organism]'
        count, ids = self.search_ids("gene", term, retmax=10)
        return GeneSearchResult(count > 0, count, ids)

    def fetch_fasta(self, ids: list[str]) -> list[tuple[str, str]]:
        if not ids:
            return []
        text = self._get("efetch.fcgi", db="nuccore", id=",".join(ids), rettype="fasta", retmode="text")
        records, header, parts = [], None, []
        for line in text.splitlines():
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(parts)))
                header, parts = line[1:].strip(), []
            elif header:
                parts.append(line.strip().upper())
        if header is not None:
            records.append((header, "".join(parts)))
        return records

    def search_and_fetch(self, term: str, retmax: int = 10) -> list[tuple[str, str]]:
        _, ids = self.search_ids("nuccore", term, retmax)
        return self.fetch_fasta(ids)
