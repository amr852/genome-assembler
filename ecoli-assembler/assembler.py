# python3
# Problem 6: Final project - de Bruijn graph assembler (baseline)
# - read-pairs are split into two single reads (pair distance is ignored)
# - k-mers packed into uint64, counted in hash partitions (fits in 512MB for E. coli X)
# - adaptive solidity threshold drops erroneous k-mers (errors -> count 1)
# - contigs = maximal non-branching paths (unitigs) of the solid k-mer graph
import sys
from collections import Counter, defaultdict
import numpy as np

K = 31
LUT = np.full(256, 4, dtype=np.uint8)
for i, ch in enumerate(b"ACGT"):
    LUT[ch] = i
BASES = "ACGT"

def read_input(B=20000):
    """Return (raw_reads_or_None, encoded_chunks, max_len).
    Reads are packed into uint8 matrices (A,C,G,T -> 0..3, padding/N -> 4) as they
    stream in, so 1.4M reads cost ~140MB instead of ~190MB of Python bytes."""
    stream = sys.stdin.buffer
    stream.readline()                         # t
    buf, chunks, raw, max_len = [], [], [], 0
    def flush():
        L = max(len(r) for r in buf)
        M = np.full((len(buf), L), 4, dtype=np.uint8)
        for i, r in enumerate(buf):
            M[i, :len(r)] = LUT[np.frombuffer(r, dtype=np.uint8)]
        chunks.append(M)
        buf.clear()
    for ln in stream:
        ln = ln.strip()
        if not ln:
            continue
        parts = ln.split(b"|")[:2] if b"|" in ln else [ln]
        for r in parts:
            if len(raw) < 50:
                raw.append(r)
            max_len = max(max_len, len(r))
            buf.append(r)
            if len(buf) >= B:
                flush()
    if buf:
        flush()
    n = sum(len(c) for c in chunks)
    return (raw if n < 50 else None), chunks, max_len

# ---------- tiny-input fallback (e.g. the sample): reads themselves are edges ----------
def small_assembly(reads):
    reads = [r.decode() for r in reads]
    k = min(len(r) for r in reads)
    edges = []
    for r in reads:
        for i in range(len(r) - k + 1):
            edges.append(r[i:i + k])
    out, indeg, outdeg = defaultdict(list), Counter(), Counter()
    for e in edges:
        a, b = e[:-1], e[1:]
        out[a].append(b); outdeg[a] += 1; indeg[b] += 1
    nodes = set(indeg) | set(outdeg)
    one = lambda v: indeg[v] == 1 and outdeg[v] == 1
    used = Counter()
    contigs = []
    for v in nodes:
        if not one(v):
            for w in out[v]:
                path = v + w[-1]
                used[(v, w)] += 1
                while one(w):
                    x = out[w][0]
                    used[(w, x)] += 1
                    path += x[-1]; w = x
                contigs.append(path)
    # isolated cycles
    for v in nodes:
        if one(v) and used[(v, out[v][0])] == 0:
            path, w = v, v
            while used[(w, out[w][0])] == 0:
                x = out[w][0]; used[(w, x)] += 1; path += x[-1]; w = x
            contigs.append(path)
    return sorted(contigs)

# ---------- main path ----------
def kmer_codes(M, k):
    L = M.shape[1]
    if L < k:
        return np.empty(0, dtype=np.uint64)
    W = L - k + 1
    bad = np.zeros((M.shape[0], L + 1), dtype=np.int16)
    np.cumsum(M == 4, axis=1, out=bad[:, 1:])
    valid = (bad[:, k:] - bad[:, :W]) == 0
    del bad
    acc = np.zeros((M.shape[0], W), dtype=np.uint64)
    for i in range(k):
        acc <<= np.uint64(2)
        acc |= (M[:, i:i + W] & 3).astype(np.uint64)
    return acc[valid]

def part_of(codes, P):
    return ((codes * np.uint64(0x9E3779B97F4A7C15)) >> np.uint64(58)) % np.uint64(P)

def pick_threshold(cnt):
    # k-mer spectrum: erroneous k-mers pile up at low counts, genuine ones form a
    # peak near the k-mer coverage. Cut at the first valley of the histogram.
    if len(cnt) == 0:
        return 1
    h = np.bincount(np.minimum(cnt, 1000))
    if h[1] < 0.15 * h[1:].sum():          # essentially error-free
        return 1
    for c in range(2, min(len(h) - 1, 60)):
        if h[c] <= h[c + 1]:
            return c
    return 2

def solid_kmers(batches, k, P=24):
    # gather per-partition codes in one encoding pass per partition (memory bounded)
    result, threshold = [], None
    for p in range(P):
        chunks = []
        for b in batches:
            c = kmer_codes(b, k)
            chunks.append(c[part_of(c, P) == p])
        allc = np.concatenate(chunks) if chunks else np.empty(0, np.uint64)
        del chunks
        u, cnt = np.unique(allc, return_counts=True)
        del allc
        if threshold is None:
            threshold = pick_threshold(cnt)
            sys.stderr.write("solidity threshold %d\n" % threshold)
        result.append(u[cnt >= threshold])
    return np.sort(np.concatenate(result))

def decode(code, k):
    s = []
    for _ in range(k):
        s.append(BASES[code & 3]); code >>= 2
    return "".join(reversed(s))

def unitigs(solid, k):
    n = len(solid)
    mask = np.uint64((1 << (2 * (k - 1))) - 1)
    pref = solid >> np.uint64(2)            # sorted, because solid is sorted
    suf = solid & mask
    suf_sorted = np.sort(suf)
    def outdeg(nodes):
        return np.searchsorted(pref, nodes, "right") - np.searchsorted(pref, nodes, "left")
    def indeg(nodes):
        return np.searchsorted(suf_sorted, nodes, "right") - np.searchsorted(suf_sorted, nodes, "left")
    simple_suf = (indeg(suf) == 1) & (outdeg(suf) == 1)
    simple_pref = (indeg(pref) == 1) & (outdeg(pref) == 1)
    nxt = np.where(simple_suf, np.searchsorted(pref, suf, "left"), -1).astype(np.int64)
    del pref, suf, suf_sorted, simple_suf
    starts = np.nonzero(~simple_pref)[0]
    del simple_pref
    codes = solid
    seen = bytearray(n)
    contigs = []
    def walk(i):
        seq = [decode(int(codes[i]), k)]
        seen[i] = 1
        j = int(nxt[i])
        while j != -1 and not seen[j]:
            seen[j] = 1
            seq.append(BASES[int(codes[j]) & 3])
            j = int(nxt[j])
        return "".join(seq)
    for i in starts.tolist():
        if not seen[i]:
            contigs.append(walk(i))
    for i in range(n):                       # leftover isolated cycles
        if not seen[i]:
            contigs.append(walk(i))
    return contigs

def main():
    raw, chunks, max_len = read_input()
    if not chunks:
        return
    if raw is not None or max_len < K:
        contigs = small_assembly(raw if raw is not None else [])
    else:
        solid = solid_kmers(chunks, K)
        del chunks
        contigs = unitigs(solid, K)
        del solid
        contigs = [c for c in contigs if len(c) >= 2 * K] or contigs   # drop tip/bubble scraps
        contigs.sort(key=len, reverse=True)
    out = sys.stdout
    for i, c in enumerate(contigs, 1):
        out.write(">CONTIG%d\n%s\n" % (i, c))

if __name__ == "__main__":
    main()
