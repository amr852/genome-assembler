# python3
# Problem 4: Tip Removal (k = 15)
# Iteratively strip edges leaving in-degree-0 nodes and entering out-degree-0 nodes;
# count removed edges. Removing one tip can expose the next, so use a work queue.
import sys
from collections import deque

K = 15

def main():
    reads = [l.strip() for l in sys.stdin if l.strip()]
    kmers = set()
    for r in reads:
        for i in range(len(r) - K + 1):
            kmers.add(r[i:i + K])
    ids = {}
    def nid(s):
        if s not in ids:
            ids[s] = len(ids)
        return ids[s]
    out_adj, in_adj = {}, {}
    for km in kmers:
        a, b = nid(km[:-1]), nid(km[1:])
        out_adj.setdefault(a, set()).add(b)
        in_adj.setdefault(b, set()).add(a)
    n = len(ids)
    for v in range(n):
        out_adj.setdefault(v, set()); in_adj.setdefault(v, set())

    removed = 0
    q = deque(range(n))
    while q:
        v = q.popleft()
        if not in_adj[v] and out_adj[v]:
            for w in list(out_adj[v]):
                out_adj[v].discard(w); in_adj[w].discard(v); removed += 1
                q.append(w)
        elif not out_adj[v] and in_adj[v]:
            for u in list(in_adj[v]):
                in_adj[v].discard(u); out_adj[u].discard(v); removed += 1
                q.append(u)
    print(removed)

if __name__ == "__main__":
    main()
