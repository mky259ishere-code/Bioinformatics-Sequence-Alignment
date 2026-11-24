"""Minimal NCBI BLASTN URL API client."""
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass


@dataclass(frozen=True)
class BlastHit:
    accession: str
    title: str
    identity: float
    coverage: float
    evalue: str
    alignment_length: int
    link: str


class BlastClient:
    URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

    def _request(self, params: dict[str, str]) -> str:
        encoded = urllib.parse.urlencode(params).encode()
        request = urllib.request.Request(
            self.URL,
            data=encoded,
            headers={"User-Agent": "BioinformaticsSequenceAlignment/0.2"},
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            return response.read().decode("utf-8", errors="replace")

    def submit(self, sequence: str, database: str = "nt", max_hits: int = 20) -> str:
        response = self._request({
            "CMD": "Put",
            "PROGRAM": "blastn",
            "DATABASE": database,
            "MEGABLAST": "on",
            "QUERY": f">user_query\n{sequence}",
            "HITLIST_SIZE": str(max_hits),
            "FORMAT_TYPE": "XML2",
        })
        rid = next((line.split("=", 1)[1].strip() for line in response.splitlines() if line.startswith(" RID =")), None)
        if not rid:
            raise RuntimeError("NCBI did not return a BLAST request ID.")
        return rid

    def fetch_xml(self, rid: str, timeout: int = 600) -> str:
        started = time.monotonic()
        while time.monotonic() - started < timeout:
            response = self._request({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML2"})
            if "Status=WAITING" in response or "Status=UNKNOWN" in response:
                time.sleep(20)
                continue
            if "Status=FAILED" in response:
                raise RuntimeError("NCBI BLAST reported that the search failed.")
            if "<BlastOutput" in response or "<BlastOutput2" in response:
                return response
            time.sleep(20)
        raise TimeoutError("NCBI BLAST timed out. Try again later.")

    def search(self, sequence: str, max_hits: int = 20) -> list[BlastHit]:
        rid = self.submit(sequence, max_hits=max_hits)
        xml = self.fetch_xml(rid)
        root = ET.fromstring(xml)
        hits = []
        for hit in root.findall(".//Hit"):
            accession = hit.findtext("Hit_accession") or hit.findtext("Hit_id") or "unknown"
            title = hit.findtext("Hit_def") or accession
            hsp = hit.find(".//Hsp")
            if hsp is None:
                continue
            length = int(hsp.findtext("Hsp_align-len", "0"))
            identities = int(hsp.findtext("Hsp_identity", "0"))
            query_from = int(hsp.findtext("Hsp_query-from", "1"))
            query_to = int(hsp.findtext("Hsp_query-to", "1"))
            coverage = 100.0 * (abs(query_to - query_from) + 1) / len(sequence)
            identity = 100.0 * identities / length if length else 0.0
            hits.append(BlastHit(
                accession, title, identity, coverage,
                hsp.findtext("Hsp_evalue", "unknown"), length,
                f"https://www.ncbi.nlm.nih.gov/nuccore/{urllib.parse.quote(accession)}",
            ))
        return hits
