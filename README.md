# Genome Assembler

Rebuilding genomes from millions of short, error-prone DNA reads, using graph algorithms written from scratch in Python.

This repo has two complete assembly pipelines, **overlap graphs** and **de Bruijn graphs**, plus the graph algorithms behind them (Eulerian cycles, max-flow circulations, tip and bubble detection). It grows from a 5 kb viral genome (φX174) up to a memory-bounded assembler for **~1.4M reads of an E. coli–scale genome under 512 MB of RAM**.

---

## Highlights

| | |
|---|---|
| **Exact reconstruction** | Both φX174 assemblers rebuild the full 5,386 bp circular genome exactly, even with a 1% sequencing error rate |
| **Memory-bounded k-mer counting** | 31-mers packed into `uint64` (2 bits/base) and counted in 24 hash partitions with NumPy, so a 1.4M-read dataset fits in 512 MB |
| **Error correction without a reference** | Picks the cutoff for erroneous k-mers automatically from the valley in the k-mer count histogram |
| **Graph cleanup** | Iterative tip removal and bubble detection, then unitig (non-branching path) construction to build contigs |
| **Tested at scale** | Simulated-read stress harness that checks the output is an exact rotation of the true circular genome |

## Benchmarks

Simulated datasets, random genomes, 100 bp reads, single-core Python 3:

| Assembler | Input | Result | Time |
|---|---|---|---|
| Overlap graph, error-free | 1,618 reads, 5,386 bp genome | **exact genome**, 5/5 trials | < 0.1 s |
| Overlap graph, error-prone (1 mismatch per read) | 1,618 reads | **exact genome**, 5/5 trials | 0.3 s |
| de Bruijn, error-prone (k = 15) | 1,618 reads | **exact genome**, 3/3 trials | < 0.1 s |
| de Bruijn from k-mer composition | 5,386 10-mers | exact genome | < 0.1 s |
| E. coli assembler (k = 31) | 40,000 read pairs, 200 kb genome, ~40× coverage, errors | **1 contig, N50 = 199.8 kb** (99.9% of genome), error-free | 8 s |

Reproduce: `cd overlap-graph && python3 stress_test.py 10`

---

## How it works

### 1. Overlap-graph assembly (`overlap-graph/`)
Every read is a node. An edge `a → b` is weighted by the longest suffix of `a` that matches a prefix of `b`.

- **Fast overlap search:** a k-mer index (k = 12) finds candidate pairs, so the assembler never compares all pairs of reads (avoids O(n²)).
- **Error tolerance:** for error-prone reads, overlaps are found on consistent *diagonals* of shared k-mers and accepted up to a ~5% mismatch budget.
- **Path finding:** a greedy Hamiltonian path takes edges heaviest first, and a union-find structure stops it from closing a cycle too early.
- **Consensus:** each read is laid out at its offset around the circular genome, and a per-position majority vote cancels out the sequencing errors.

### 2. de Bruijn graph assembly (`de-bruijn/`)
Every k-mer is an edge `prefix(k−1) → suffix(k−1)`. A circular genome spells an **Eulerian cycle**, a walk that uses every edge exactly once.

- **Eulerian cycles** with iterative Hierholzer, O(V + E), no recursion limits.
- **k-universal strings:** a de Bruijn sequence of every binary k-mer, for k up to 14.
- **Choosing k:** the largest k whose graph is balanced and connected, so it has exactly one Eulerian cycle.
- **Tip removal:** iteratively strips dead-end branches caused by read errors, using a work queue, because removing one tip can expose the next.
- **Bubble detection:** counts pairs of short paths between the same two nodes that share no node in between (typically sequencing errors or SNPs).
- **Error-prone assembly:** drops low-coverage k-mers, removes tips, then follows the highest-coverage edges around the cycle, which smooths bubbles over.

### 3. E. coli-scale assembler (`ecoli-assembler/`)
The final project: turn millions of reads (single or paired) into FASTA contigs.

```
reads ─► 2-bit pack into uint8 matrices ─► 31-mers as uint64
      ─► hash-partitioned counting (np.unique per partition)
      ─► error cutoff from the k-mer histogram valley
      ─► solid k-mer graph ─► unitigs via sorted-array searchsorted
      ─► drop short scraps ─► FASTA contigs
```

- Reads are packed as they stream in (~140 MB for 1.4M reads instead of ~190 MB of Python byte strings).
- k-mers containing `N` are masked with a cumulative-sum window trick, so the inner loop has no Python-level branching.
- Graph degrees come from `searchsorted` on sorted prefix and suffix arrays. No hash map holding millions of entries is ever built.

### 4. Graph algorithms (`graph-algorithms/`)
- **Circulation with lower bounds:** reduced to max-flow (shift each edge's lower bound into node demands, add a super source and sink) and solved with **Dinic's algorithm**.
- **Eulerian cycle:** the general directed-graph version used by the de Bruijn tools.
- **Puzzle assembly:** a backtracking warm-up that places n×n square pieces with matching edge colors, the same "assemble fragments by overlap" idea in 2D.

---

## Repository layout

```
genome-assembler/
├── overlap-graph/
│   ├── assembler_error_free.py    # φX174, exact overlaps
│   ├── assembler_error_prone.py   # φX174, tolerant overlaps + consensus
│   └── stress_test.py             # simulated-read correctness/timing harness
├── de-bruijn/
│   ├── phi_x174_kmer.py           # genome from k-mer composition
│   ├── k_universal_string.py      # binary de Bruijn sequences
│   ├── optimal_kmer.py            # pick k giving a unique Eulerian cycle
│   ├── tip_removal.py
│   ├── bubble_detection.py
│   └── phi_x174_assembly.py       # error-prone reads → circular genome
├── ecoli-assembler/
│   └── assembler.py               # NumPy, memory-bounded, FASTA contigs
└── graph-algorithms/
    ├── eulerian_cycle.py
    ├── circulation.py             # lower-bound flows via Dinic
    └── puzzle_assembly.py
```

## Usage

Every program reads from stdin and writes to stdout.

```bash
pip install -r requirements.txt     # numpy (only needed by ecoli-assembler)

python3 overlap-graph/assembler_error_prone.py < reads.txt > genome.txt
python3 de-bruijn/phi_x174_assembly.py        < reads.txt > genome.txt
python3 ecoli-assembler/assembler.py          < reads.txt > contigs.fasta
echo 4 | python3 de-bruijn/k_universal_string.py      # → 0000110010111101 (one valid answer)
```

**Input formats**
- Assemblers: one read per line. The E. coli assembler's first line is the read count, and paired reads look like `READ1|READ2|d`.
- `bubble_detection.py`: first line `k t`, then reads.
- `circulation.py`: `n m`, then `m` lines of `u v lower capacity`.

## Limitations and next steps

- **Pair distances are ignored:** paired reads are split into two single reads. Scaffolding with the insert size would join contigs across repeats.
- **Bubbles are smoothed over, not resolved:** the φX174 de Bruijn walk follows the higher-coverage branch, and the E. coli assembler drops short contigs instead of merging the two branches of a bubble.
- **Only one DNA strand is used:** the simulated data comes from a single strand. Real data needs canonical (strand-independent) k-mers.
- Possible upgrades: tune k automatically per dataset, write the hot loops in Rust or C via PyO3, and compare NGA50 against SPAdes.

## Background

Built while working through the *Genome Assembly Programming Challenge* (UC San Diego, Coursera). The algorithms follow Compeau & Pevzner, *Bioinformatics Algorithms*.

## License

MIT
