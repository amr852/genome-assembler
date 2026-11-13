# python3
"""Problem 4: Assemble the circular phi X174 genome from its k-mer composition.

Each k-mer is an edge prefix(k-1) -> suffix(k-1) in the de Bruijn graph.
An Eulerian cycle through all edges spells the circular genome; its length
equals the number of k-mers.
"""
import sys


def assemble(kmers):
    ids = {}
    labels = []

    def node(s):
        if s not in ids:
            ids[s] = len(labels)
            labels.append(s)
        return ids[s]

    adj = []
    for kmer in kmers:
        u, v = node(kmer[:-1]), node(kmer[1:])
        while len(adj) < len(labels):
            adj.append([])
        adj[u].append(v)
    while len(adj) < len(labels):
        adj.append([])

    ptr = [0] * len(labels)
    stack = [0]
    path = []
    while stack:
        v = stack[-1]
        if ptr[v] < len(adj[v]):
            stack.append(adj[v][ptr[v]])
            ptr[v] += 1
        else:
            path.append(stack.pop())
    path.reverse()               # len(kmers) + 1 nodes, first == last
    # Circular genome: first character of each node along the cycle.
    return "".join(labels[v][0] for v in path[:-1])


def main():
    kmers = [l.strip() for l in sys.stdin if l.strip()]
    print(assemble(kmers))


if __name__ == "__main__":
    main()
