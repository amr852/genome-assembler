# python3
# Problem 2: Selecting the Optimal k-mer Size
# Largest k whose de Bruijn graph (edges = distinct k-mers, nodes = (k-1)-mers)
# is balanced (indeg == outdeg everywhere) and connected -> a single Eulerian cycle.
import sys
from collections import defaultdict

def is_eulerian(reads, k):
    kmers = set()
    for r in reads:
        for i in range(len(r) - k + 1):
            kmers.add(r[i:i + k])
    indeg, outdeg = defaultdict(int), defaultdict(int)
    adj = defaultdict(list)
    for km in kmers:
        a, b = km[:-1], km[1:]
        outdeg[a] += 1; indeg[b] += 1
        adj[a].append(b)
    nodes = set(indeg) | set(outdeg)
    for v in nodes:
        if indeg[v] != outdeg[v]:
            return False
    # balanced + weakly connected => strongly connected => Eulerian cycle
    und = defaultdict(list)
    for a in adj:
        for b in adj[a]:
            und[a].append(b); und[b].append(a)
    start = next(iter(nodes))
    seen = {start}; stack = [start]
    while stack:
        u = stack.pop()
        for w in und[u]:
            if w not in seen:
                seen.add(w); stack.append(w)
    return len(seen) == len(nodes)

def main():
    reads = [l.strip() for l in sys.stdin if l.strip()]
    L = min(len(r) for r in reads)
    for k in range(L, 1, -1):
        if is_eulerian(reads, k):
            print(k)
            return
    print(2)

if __name__ == "__main__":
    main()
