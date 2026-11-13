# python3
"""Problem 2: Eulerian cycle in a directed graph — iterative Hierholzer, O(n + m)."""
import sys


def eulerian_cycle(n, edges):
    adj = [[] for _ in range(n + 1)]
    indeg = [0] * (n + 1)
    for u, v in edges:
        adj[u].append(v)
        indeg[v] += 1

    # Graph is guaranteed strongly connected, so balance is the only check.
    for v in range(1, n + 1):
        if len(adj[v]) != indeg[v]:
            return None

    ptr = [0] * (n + 1)          # next unused edge per vertex
    stack = [edges[0][0]]
    path = []
    while stack:
        v = stack[-1]
        if ptr[v] < len(adj[v]):
            stack.append(adj[v][ptr[v]])
            ptr[v] += 1
        else:
            path.append(stack.pop())
    path.reverse()               # m + 1 vertices, first == last
    return path[:-1]


def main():
    data = sys.stdin.buffer.read().split()
    n, m = int(data[0]), int(data[1])
    nums = list(map(int, data[2:2 + 2 * m]))
    edges = list(zip(nums[0::2], nums[1::2]))
    cycle = eulerian_cycle(n, edges)
    if cycle is None:
        print(0)
    else:
        print(1)
        print(" ".join(map(str, cycle)))


if __name__ == "__main__":
    main()
