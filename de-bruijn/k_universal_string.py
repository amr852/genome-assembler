# python3
"""Problem 3: k-universal circular binary string.

De Bruijn graph: nodes = all (k-1)-bit strings (as ints), each node v has
edges to (v << 1 | b) & mask for b in {0,1}. Every edge is a k-mer, so an
Eulerian cycle spells a string containing each k-mer exactly once.
"""
import sys


def k_universal(k):
    if k == 1:
        return "01"
    nodes = 1 << (k - 1)
    mask = nodes - 1
    ptr = [0] * nodes            # next unused out-bit (0, then 1) per node
    stack = [0]
    path = []
    while stack:
        v = stack[-1]
        if ptr[v] < 2:
            b = ptr[v]
            ptr[v] += 1
            stack.append(((v << 1) | b) & mask)
        else:
            path.append(stack.pop())
    path.reverse()               # 2^k + 1 nodes, first == last
    # Each step appends the new low bit; take the last bit of every node after the first.
    return "".join(str(v & 1) for v in path[1:])


def main():
    k = int(sys.stdin.readline())
    print(k_universal(k))


if __name__ == "__main__":
    main()
