import tempfile
import unittest
from pathlib import Path

from src.alignment import smith_waterman
from src.fasta import read_fasta


class AlignmentTests(unittest.TestCase):
    def test_local_alignment_finds_shared_region(self):
        result = smith_waterman("GGTACGTTAGC", "TTACGTTAGCT")
        self.assertGreaterEqual(result.score, 16)
        self.assertGreaterEqual(result.identity_percent, 90.0)

    def test_fasta_parser_reads_multiple_records(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.fasta"
            path.write_text(">one\\nACGT\\n>two\\nTTAA\\n", encoding="utf-8")
            self.assertEqual(read_fasta(path), [("one", "ACGT"), ("two", "TTAA")])


if __name__ == "__main__":
    unittest.main()
