# python3
# Problem 1: Finding a Circulation in a Network (lower bounds -> max flow)
import sys
from collections import deque

def main():
    data = sys.stdin.read().split()
    n, m = int(data[0]), int(data[1])
    S, T = 0, n + 1
    N = n + 2
    # adjacency with edge arrays: to, cap, rev index via pairing (e ^ 1)
    to, cap, head, nxt = [], [], [-1] * N, []

    def add(u, v, c):
        to.append(v); cap.append(c); nxt.append(head[u]); head[u] = len(to) - 1
        to.append(u); cap.append(0); nxt.append(head[v]); head[v] = len(to) - 1
        return len(to) - 2

    excess = [0] * N
    edge_ids, lows = [], []
    idx = 2
    for _ in range(m):
        u, v, l, c = (int(x) for x in data[idx:idx + 4]); idx += 4
        edge_ids.append(add(u, v, c - l))
        lows.append(l)
        excess[v] += l
        excess[u] -= l
    need = 0
    for v in range(1, n + 1):
        if excess[v] > 0:
            add(S, v, excess[v]); need += excess[v]
        elif excess[v] < 0:
            add(v, T, -excess[v])

    # Dinic
    flow = 0
    while True:
        level = [-1] * N
        level[S] = 0
        q = deque([S])
        while q:
            u = q.popleft()
            e = head[u]
            while e != -1:
                if cap[e] > 0 and level[to[e]] < 0:
                    level[to[e]] = level[u] + 1
                    q.append(to[e])
                e = nxt[e]
        if level[T] < 0:
            break
        it = head[:]

        def dfs(u, f):
            if u == T:
                return f
            while it[u] != -1:
                e = it[u]
                v = to[e]
                if cap[e] > 0 and level[v] == level[u] + 1:
                    d = dfs(v, min(f, cap[e]))
                    if d:
                        cap[e] -= d; cap[e ^ 1] += d
                        return d
                it[u] = nxt[e]
            return 0

        while True:
            f = dfs(S, float('inf'))
            if not f:
                break
            flow += f

    if flow != need:
        print("NO")
        return
    out = ["YES"]
    for e, l in zip(edge_ids, lows):
        out.append(str(l + cap[e ^ 1]))  # flow on edge = reverse residual
    print("\n".join(out))

if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
