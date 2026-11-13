# python3
"""Stress test: simulate phi X174-like datasets and check the assemblers.
Usage: python3 stress_test.py [trials]"""
import random, subprocess, sys, time

GENOME_LEN, N_READS, READ_LEN = 5386, 1618, 100


def make_genome(rng):
    return "".join(rng.choice("ACGT") for _ in range(GENOME_LEN))


def make_reads(genome, rng, errors):
    g2 = genome + genome[:READ_LEN]
    reads = []
    for _ in range(N_READS):
        s = rng.randrange(GENOME_LEN)
        r = list(g2[s:s + READ_LEN])
        if errors:
            i = rng.randrange(READ_LEN)
            r[i] = rng.choice([c for c in "ACGT" if c != r[i]])
        reads.append("".join(r))
    return sorted(reads)


def is_rotation(a, b):
    return len(a) == len(b) and a in b + b


def hamming_best_rotation(out, genome):
    if len(out) != len(genome):
        return None
    doubled = out + out
    # align via a 20-mer anchor
    anchor = genome[:20]
    best = None
    pos = doubled.find(anchor)
    while pos != -1 and pos < len(out):
        rot = doubled[pos:pos + len(out)]
        d = sum(x != y for x, y in zip(rot, genome))
        best = d if best is None else min(best, d)
        pos = doubled.find(anchor, pos + 1)
    return best


def run(script, data):
    t = time.time()
    out = subprocess.run([sys.executable, script], input=data, capture_output=True, text=True).stdout.strip()
    return out, time.time() - t


trials = int(sys.argv[1]) if len(sys.argv) > 1 else 10
ok = {1: 0, 2: 0}
tmax = {1: 0.0, 2: 0.0}
for seed in range(trials):
    rng = random.Random(seed)
    g = make_genome(rng)
    for prob, script in ((1, "assembler_error_free.py"), (2, "assembler_error_prone.py")):
        reads = make_reads(g, rng, errors=(prob == 2))
        out, dt = run(script, "\n".join(reads) + "\n")
        tmax[prob] = max(tmax[prob], dt)
        good = is_rotation(out, g)
        ok[prob] += good
        if not good:
            print("seed %d P%d FAIL: len %d vs %d, hamming %s" % (seed, prob, len(out), len(g), hamming_best_rotation(out, g)))
for p in (1, 2):
    print("P%d: %d/%d exact rotation, max time %.2fs" % (p, ok[p], trials, tmax[p]))
