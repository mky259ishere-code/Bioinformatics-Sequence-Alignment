"""Asynchronous NCBI BLASTN URL API client."""
import re
import ssl
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


def system_ssl_context() -> ssl.SSLContext:
    try:
        import truststore
        return truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    except ImportError:
        return ssl.create_default_context()


class BlastClient:
    URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    RESULT_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID={}"

    def _request(self, params: dict[str, str]) -> str:
        request = urllib.request.Request(self.URL, data=urllib.parse.urlencode(params).encode(),
                                         headers={"User-Agent": "BioinformaticsSequenceAlignment/0.3"})
        with urllib.request.urlopen(request, timeout=60, context=system_ssl_context()) as response:
            return response.read().decode("utf-8", errors="replace")

    def submit(self, sequence: str, max_hits: int = 20) -> str:
        response = self._request({
            "CMD": "Put", "PROGRAM": "blastn", "DATABASE": "nt",
            "MEGABLAST": "on", "SHORT_QUERY_ADJUST": "true",
            "QUERY": f">user_query\n{sequence}",
            "HITLIST_SIZE": str(max_hits), "FORMAT_TYPE": "XML2",
        })
        match = re.search(r"\bRID\s*=\s*([A-Za-z0-9_-]+)", response)
        if match:
            return match.group(1)
        raise RuntimeError("NCBI did not return a BLAST request ID. Response: " + response[:300].replace("\n", " "))

    def fetch_xml(self, rid: str) -> str:
        response = self._request({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML2"})
        if "Status=WAITING" in response or "Status=UNKNOWN" in response:
            raise RuntimeError("Results are not ready yet.")
        if "Status=FAILED" in response:
            raise RuntimeError("NCBI BLAST reported that the search failed.")
        if "<BlastOutput" not in response and "<BlastOutput2" not in response:
            raise RuntimeError("Results are not ready yet.")
        return response

    def parse_results(self, xml: str, query_length: int) -> list[BlastHit]:
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
            q_from = int(hsp.findtext("Hsp_query-from", "1"))
            q_to = int(hsp.findtext("Hsp_query-to", "1"))
            coverage = 100.0 * (abs(q_to - q_from) + 1) / query_length
            identity = 100.0 * identities / length if length else 0.0
            hits.append(BlastHit(accession, title, identity, coverage,
                                 hsp.findtext("Hsp_evalue", "unknown"), length,
                                 f"https://www.ncbi.nlm.nih.gov/nuccore/{urllib.parse.quote(accession)}"))
        return hits
