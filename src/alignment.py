"""Smith-Waterman local alignment for DNA sequences."""

from dataclasses import dataclass


@dataclass(frozen=True)
class AlignmentResult:
    score: int
    aligned_query: str
    aligned_reference: str
    query_start: int
    query_end: int
    reference_start: int
    reference_end: int

    @property
    def alignment_length(self) -> int:
        return len(self.aligned_query)

    @property
    def matches(self) -> int:
        return sum(a == b and a != "-" for a, b in zip(self.aligned_query, self.aligned_reference))

    @property
    def identity_percent(self) -> float:
        return 100.0 * self.matches / self.alignment_length if self.alignment_length else 0.0


def smith_waterman(query: str, reference: str, match: int = 2, mismatch: int = -1, gap: int = -2) -> AlignmentResult:
    """Find the highest-scoring local alignment between two DNA sequences."""
    query, reference = query.upper(), reference.upper()
    rows, columns = len(query) + 1, len(reference) + 1
    scores = [[0] * columns for _ in range(rows)]
    best_score = 0
    best_cell = (0, 0)

    for i in range(1, rows):
        for j in range(1, columns):
            diagonal = scores[i - 1][j - 1] + (match if query[i - 1] == reference[j - 1] else mismatch)
            upward = scores[i - 1][j] + gap
            left = scores[i][j - 1] + gap
            scores[i][j] = max(0, diagonal, upward, left)
            if scores[i][j] > best_score:
                best_score, best_cell = scores[i][j], (i, j)

    i, j = best_cell
    query_end, reference_end = i, j
    aligned_query: list[str] = []
    aligned_reference: list[str] = []

    while i > 0 and j > 0 and scores[i][j] > 0:
        current = scores[i][j]
        if current == scores[i - 1][j - 1] + (match if query[i - 1] == reference[j - 1] else mismatch):
            aligned_query.append(query[i - 1]); aligned_reference.append(reference[j - 1]); i -= 1; j -= 1
        elif current == scores[i - 1][j] + gap:
            aligned_query.append(query[i - 1]); aligned_reference.append("-"); i -= 1
        else:
            aligned_query.append("-"); aligned_reference.append(reference[j - 1]); j -= 1

    return AlignmentResult(best_score, "".join(reversed(aligned_query)), "".join(reversed(aligned_reference)), i, query_end, j, reference_end)
