"""BLAST-like seed-and-extend search for nucleotide sequences."""
from dataclasses import dataclass
from .alignment import smith_waterman


@dataclass(frozen=True)
class SearchHit:
    name: str
    score: int
    identity: float
    seeds: int
    span: int
    aligned_query: str
    aligned_reference: str


def _seed_positions(sequence: str, k: int) -> dict[str, list[int]]:
    return {sequence[i:i+k]: [i] for i in range(len(sequence)-k+1)}


def search_references(query: str, references: list[tuple[str, str]], k: int = 11, min_seeds: int = 1) -> list[SearchHit]:
    if k < 2 or k > len(query):
        raise ValueError("k-mer length must be between 2 and the query length")
    query = query.upper()
    query_seeds = _seed_positions(query, k)
    hits: list[SearchHit] = []
    for name, reference in references:
        index = {}
        for i in range(len(reference)-k+1):
            index.setdefault(reference[i:i+k], []).append(i)
        seed_count = sum(len(index.get(seed, [])) for seed in query_seeds)
        if seed_count < min_seeds:
            continue
        result = smith_waterman(query, reference)
        hits.append(SearchHit(name, result.score, result.identity_percent, seed_count,
                              result.alignment_length, result.aligned_query, result.aligned_reference))
    return sorted(hits, key=lambda h: (h.score, h.identity, h.seeds), reverse=True)
