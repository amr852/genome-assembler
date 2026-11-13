# python3
"""
Programming Assignment 1, Problem 1:
Assembling phi X174 from error-free reads using an overlap graph.

Approach
--------
1. Deduplicate reads (identical reads add nothing to an error-free assembly).
2. Index the k-mer prefix of every read (k = 12). An overlap of length L from
   read a into read b means b's prefix k-mer appears in a at position len(a) - L,
   so scanning a's k-mers left-to-right yields candidates in decreasing L and we
   only verify pairs that actually share a k-mer (no all-pairs comparison).
3. Overlap graph: edge a -> b weighted by the longest exact suffix/prefix
   overlap. Keep the TOP_EDGES heaviest out-edges per read.
4. Greedy Hamiltonian path over those edges, spell it, then trim the
   last -> first overlap because the genome is circular.
"""
import sys

K = 12
TOP_EDGES = 5


def read_input():
    seen, uniq = set(), []
    for line in sys.stdin:
        r = line.strip()
        if r and r not in seen:
            seen.add(r)
            uniq.append(r)
    return uniq


def build_prefix_index(reads, k):
    index = {}
    for idx, r in enumerate(reads):
        index.setdefault(r[:k], []).append(idx)
    return index


def out_edges(a_idx, reads, index, k, limit):
    """Up to `limit` heaviest exact overlaps a -> b as (L, b)."""
    a = reads[a_idx]
    n = len(a)
    res, taken = [], set()
    for i in range(1, n - k + 1):
        cands = index.get(a[i:i + k])
        if not cands:
            continue
        suffix = a[i:]
        for b_idx in cands:
            if b_idx != a_idx and b_idx not in taken and reads[b_idx].startswith(suffix):
                res.append((n - i, b_idx))
                taken.add(b_idx)
                if len(res) >= limit:
                    return res
    return res


def overlap_len(a, b):
    for L in range(min(len(a), len(b)) - 1, 0, -1):
        if a[-L:] == b[:L]:
            return L
    return 0

def greedy_hamiltonian_path(n, edges):
    """Global greedy: take edges by decreasing weight; each read gets at most one
    successor and one predecessor, and no edge may close a cycle early
    (union-find). Returns the ordered list of reads and the edge payloads used.
    edges: list of (weight, a, b, payload)."""
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    succ = [None] * n
    pred = [None] * n
    payload = {}
    used = 0
    edges.sort(key=lambda e: (-e[0], e[1], e[2]))
    for w, a, b, pl in edges:
        if used == n - 1:
            break
        if succ[a] is not None or pred[b] is not None:
            continue
        ra, rb = find(a), find(b)
        if ra == rb:
            continue
        parent[ra] = rb
        succ[a] = b
        pred[b] = a
        payload[(a, b)] = pl
        used += 1

    # Walk the fragments (normally exactly one) and chain them together.
    path, links = [], []
    for start in range(n):
        if pred[start] is not None:
            continue
        cur = start
        if path:
            links.append((path[-1], cur, None))   # fragment junction, no overlap known
        path.append(cur)
        while succ[cur] is not None:
            nxt = succ[cur]
            links.append((cur, nxt, payload[(cur, nxt)]))
            path.append(nxt)
            cur = nxt
    return path, links


def assemble(reads, k=K):
    if not reads:
        return ""
    if len(reads) == 1:
        return reads[0]
    k = min(k, max(1, min(len(r) for r in reads) // 2))
    index = build_prefix_index(reads, k)

    edges = []
    for a in range(len(reads)):
        for L, b in out_edges(a, reads, index, k, TOP_EDGES):
            edges.append((L, a, b, L))

    path, links = greedy_hamiltonian_path(len(reads), edges)

    pieces = [reads[path[0]]]
    for a, b, L in links:
        if L is None:
            L = overlap_len(reads[a], reads[b])
        pieces.append(reads[b][L:])
    genome = "".join(pieces)

    # Circular genome: the last read wraps around into the first read.
    wrap = overlap_len(reads[path[-1]], reads[path[0]])
    if wrap:
        genome = genome[:-wrap]
    return genome


def main():
    sys.stdout.write(assemble(read_input()) + "\n")


if __name__ == "__main__":
    main()
