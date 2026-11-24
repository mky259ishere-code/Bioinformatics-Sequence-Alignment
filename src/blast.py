"""Minimal NCBI BLASTN URL API client."""
import re
import ssl
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
                                         headers={"User-Agent": "BioinformaticsSequenceAlignment/0.2"})
        try:
            with urllib.request.urlopen(request, timeout=60, context=system_ssl_context()) as response:
                return response.read().decode("utf-8", errors="replace")
        except ssl.SSLCertVerificationError as error:
            raise RuntimeError("SSL verification failed. Run 'python -m pip install truststore'.") from error

    def submit(self, sequence: str, database: str = "nt", max_hits: int = 20) -> str:
        response = self._request({"CMD": "Put", "PROGRAM": "blastn", "DATABASE": database,
                                  "MEGABLAST": "on", "QUERY": f">user_query\n{sequence}",
                                  "HITLIST_SIZE": str(max_hits), "FORMAT_TYPE": "XML2"})
        match = re.search(r"\bRID\s*=\s*([A-Za-z0-9_-]+)", response)
        if match:
            return match.group(1)
        detail = " ".join(line.strip() for line in response.splitlines()
                          if "ERROR" in line.upper() or "MESSAGE" in line.upper())[:500]
        raise RuntimeError("NCBI did not return a BLAST request ID." + (f" Server response: {detail}" if detail else ""))

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
        raise TimeoutError(
            f"NCBI is still processing the search. Open this result link later: {self.RESULT_URL.format(rid)}"
        )

    def search(self, sequence: str, max_hits: int = 20) -> list[BlastHit]:
        root = ET.fromstring(self.fetch_xml(self.submit(sequence, max_hits=max_hits)))
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
            hits.append(BlastHit(accession, title, identity, coverage,
                                 hsp.findtext("Hsp_evalue", "unknown"), length,
                                 f"https://www.ncbi.nlm.nih.gov/nuccore/{urllib.parse.quote(accession)}"))
        return hits
