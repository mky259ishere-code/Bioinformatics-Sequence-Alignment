"""Small, rate-limited NCBI Entrez client using ESearch and EFetch."""
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


class NCBIClient:
    BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    def __init__(self, email: str, api_key: str | None = None, tool: str = "BioinformaticsSequenceAlignment"):
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

    def search_ids(self, term: str, retmax: int = 10) -> list[str]:
        xml = self._get("esearch.fcgi", db="nuccore", term=term, retmode="xml", retmax=retmax)
        return [node.text for node in ET.fromstring(xml).findall(".//Id") if node.text]

    def fetch_fasta(self, ids: list[str]) -> list[tuple[str, str]]:
        if not ids:
            return []
        text = self._get("efetch.fcgi", db="nuccore", id=",".join(ids), rettype="fasta", retmode="text")
        records = []
        header, parts = None, []
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
        return self.fetch_fasta(self.search_ids(term, retmax))
