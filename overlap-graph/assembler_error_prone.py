# python3
"""
Programming Assignment 1, Problem 2:
Assembling phi X174 from error-prone reads (1 mismatch per 100-mer) using an
overlap graph.

Approach
--------
1. Index every k-mer of every read (k = 12). With only one error per read,
   most k-mers survive, so true overlaps share k-mers on a consistent
   diagonal d = p - q (a's k-mer at p, b's at q), i.e. overlap L = len(a) - d.
2. For each read, try candidate diagonals in decreasing L and accept an overlap
   if its mismatch count is within 5% of L (+1 slack on long overlaps, since
   both reads may carry their error inside the overlap). d = 0 is allowed:
   two reads sampled from the same spot just stack. Keep the TOP_EDGES
   heaviest valid out-edges per read.
3. Greedy Hamiltonian path over the overlap graph.
4. Lay each read out at its offset around the circle (genome length = sum of
   offsets along the cycle, including last -> first) and take a per-position
   majority vote. The consensus is what cancels out the sequencing errors.
"""
import sys
from collections import defaultdict

K = 12
ERROR_RATE = 0.05
MIN_OVERLAP = 20
TOP_EDGES = 5


def read_input():
    return [line.strip() for line in sys.stdin if line.strip()]


def build_kmer_index(reads, k):
    index = defaultdict(list)
    for idx, r in enumerate(reads):
        for q in range(len(r) - k + 1):
            index[r[q:q + k]].append((idx, q))
    return index


def allowed_mismatches(L):
    return int(ERROR_RATE * L) + (1 if L >= 20 else 0)


def mismatches(a, b, d, limit):
    mm = 0
    for x, y in zip(a[d:], b):
        if x != y:
            mm += 1
            if mm > limit:
                break
    return mm


def out_edges(a_idx, reads, index, k, limit, only=None):
    """Up to `limit` heaviest tolerant overlaps a -> b as (L, b, d)."""
    a = reads[a_idx]
    n = len(a)
    diag = set()
    for p in range(n - k + 1):
        for b_idx, q in index.get(a[p:p + k], ()):
            if b_idx != a_idx and p >= q and (only is None or b_idx == only):
                diag.add((b_idx, p - q))
    min_overlap = min(MIN_OVERLAP, max(1, n // 3))
    cands = []
    for b_idx, d in diag:
        L = min(n - d, len(reads[b_idx]))
        if L >= min_overlap:
            cands.append((L, b_idx, d))
    cands.sort(key=lambda t: (-t[0], t[1]))
    res, taken = [], set()
    for L, b_idx, d in cands:
        if b_idx in taken:
            continue
        lim = allowed_mismatches(L)
        if mismatches(a, reads[b_idx], d, lim) <= lim:
            res.append((L, b_idx, d))
            taken.add(b_idx)
            if len(res) >= limit:
                break
    return res

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
    index = build_kmer_index(reads, k)
    n_reads = len(reads)

    edges = []
    for a in range(n_reads):
        for L, b, d in out_edges(a, reads, index, k, TOP_EDGES):
            edges.append((L, a, b, d))

    path, links = greedy_hamiltonian_path(n_reads, edges)

    shifts = [d if d is not None else len(reads[a]) for a, b, d in links]
    wrap = out_edges(path[-1], reads, index, k, 1, only=path[0])
    wrap_d = wrap[0][2] if wrap else len(reads[path[-1]])
    genome_len = sum(shifts) + wrap_d
    if genome_len <= 0:
        genome_len = max(len(r) for r in reads)

    # Majority-vote consensus around the circle.
    code = {"A": 0, "C": 1, "G": 2, "T": 3}
    counts = [[0, 0, 0, 0] for _ in range(genome_len)]
    offset = 0
    for j, r_idx in enumerate(path):
        if j:
            offset += shifts[j - 1]
        for t, ch in enumerate(reads[r_idx]):
            c = code.get(ch)
            if c is not None:
                counts[(offset + t) % genome_len][c] += 1
    out = []
    for col in counts:
        best = max(range(4), key=lambda c: col[c])
        out.append("ACGT"[best])
    return "".join(out)


def main():
    sys.stdout.write(assemble(read_input()) + "\n")


if __name__ == "__main__":
    main()
