# python3
# Problem 3: Bubble Detection
# Count pairs of short (<= t edges), simple paths v->w that share no nodes except v and w.
import sys
from collections import defaultdict

def main():
    lines = sys.stdin.read().split()
    k, t = int(lines[0]), int(lines[1])
    reads = lines[2:]
    kmers = set()
    for r in reads:
        for i in range(len(r) - k + 1):
            kmers.add(r[i:i + k])
    ids = {}
    def nid(s):
        if s not in ids:
            ids[s] = len(ids)
        return ids[s]
    edges = set()
    for km in kmers:
        edges.add((nid(km[:-1]), nid(km[1:])))
    n = len(ids)
    adj = [[] for _ in range(n)]
    indeg = [0] * n
    for a, b in edges:
        adj[a].append(b); indeg[b] += 1

    total = 0
    for v in range(n):
        if len(adj[v]) < 2:
            continue
        paths_to = defaultdict(list)  # w -> list of frozenset(intermediate nodes)
        # iterative DFS over simple paths of length <= t
        stack = [(v, [v])]
        while stack:
            u, path = stack.pop()
            if len(path) - 1 >= t:
                continue
            for w in adj[u]:
                if w in path:
                    continue
                np_ = path + [w]
                if indeg[w] >= 2:
                    paths_to[w].append(frozenset(np_[1:-1]))
                stack.append((w, np_))
        for w, plist in paths_to.items():
            for i in range(len(plist)):
                pi = plist[i]
                for j in range(i + 1, len(plist)):
                    if pi.isdisjoint(plist[j]):
                        total += 1
    print(total)

if __name__ == "__main__":
    main()
