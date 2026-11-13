# python3
# Problem 5: Assembling phi X174 from error-prone reads (k = 15)
# 1) count 15-mers, drop low-coverage (erroneous) ones
# 2) strip tips iteratively
# 3) walk the cycle greedily by highest coverage (this collapses bubbles)
import sys
from collections import Counter, deque

K = 15

def assemble(reads, min_cov):
    counts = Counter()
    for r in reads:
        for i in range(len(r) - K + 1):
            counts[r[i:i + K]] += 1
    kmers = {km: c for km, c in counts.items() if c >= min_cov}
    out_adj, in_adj = {}, {}
    for km in kmers:
        a, b = km[:-1], km[1:]
        out_adj.setdefault(a, set()).add(b); in_adj.setdefault(b, set()).add(a)
        out_adj.setdefault(b, set()); in_adj.setdefault(a, set())
    # tip removal
    q = deque(out_adj)
    while q:
        v = q.popleft()
        if not in_adj[v] and out_adj[v]:
            for w in list(out_adj[v]):
                out_adj[v].discard(w); in_adj[w].discard(v); q.append(w)
        elif not out_adj[v] and in_adj[v]:
            for u in list(in_adj[v]):
                in_adj[v].discard(u); out_adj[u].discard(v); q.append(u)
    live = [v for v in out_adj if out_adj[v]]
    if not live:
        return ""
    cov = lambda a, b: kmers.get(a + b[-1], 0)
    # start from the best-supported edge
    start = max(live, key=lambda v: max(cov(v, w) for w in out_adj[v]))
    used = set()
    path = [start]
    cur = start
    while True:
        cands = [w for w in out_adj[cur] if (cur, w) not in used]
        if not cands:
            break
        nxt = max(cands, key=lambda w: cov(cur, w))
        used.add((cur, nxt))
        cur = nxt
        if cur == start:
            break
        path.append(cur)
    # circular genome: one char per node
    return "".join(p[0] for p in path)

def greedy_circular(reads):
    """Tiny inputs (reads shorter than K): greedy shortest-superstring merge,
    then trim the wrap-around overlap of the circular genome."""
    def ov(a, b):
        for L in range(min(len(a), len(b)) - 1, 0, -1):
            if a.endswith(b[:L]):
                return L
        return 0
    frags = sorted(set(reads))
    frags = [f for f in frags if not any(f != g and f in g for g in frags)]
    while len(frags) > 1:
        best = (-1, 0, 1)
        for i, a in enumerate(frags):
            for j, b in enumerate(frags):
                if i != j:
                    L = ov(a, b)
                    if L > best[0]:
                        best = (L, i, j)
        L, i, j = best
        merged = frags[i] + frags[j][L:]
        frags = [f for t, f in enumerate(frags) if t not in (i, j)] + [merged]
    g = frags[0]
    for L in range(len(g) - 1, 0, -1):
        if g[:L] == g[-L:]:
            return g[:-L]
    return g


def main():
    global K
    reads = [l.strip() for l in sys.stdin if l.strip()]
    if len(reads[0]) < K:  # tiny sample: greedy overlap fallback
        print(greedy_circular(reads)); return
    best = ""
    for min_cov in (3, 2, 4):
        g = assemble(reads, min_cov)
        if 5000 <= len(g) <= 5800:
            best = g; break
        if len(g) > len(best):
            best = g
    print(best)

if __name__ == "__main__":
    main()
