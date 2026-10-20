# Index optimisation: design notes

Status: draft for discussion (2026-10-20). Nothing here is implemented yet.

## 1. Goal and constraints

Make the index smaller and classification faster, without losing accuracy, so the
whole workflow (build *and* classify) runs on a high-spec laptop:

| resource | budget |
|---|---|
| RAM | 16 GB total → index in RAM should be ≲ 8 GB (leave room for OS, reads, threads) |
| threads | 8 |
| disk | ~300 GB free |

The finished tool must still work for:

- **ONT** reads (main use case; Hall data minimum read length 500 bp), and
- **Illumina** paired reads (currently unsupported because results with w=41, k=19
  were poor — see §5.3).

ONT and Illumina **may use separate indexes** with different (w, k, FPR). Each index
must meet the laptop budget on its own.

**Two releases** (decided 2026-10-20):

| release | reference | must fit 16 GB laptop? | use |
|---|---|---|---|
| **small** | `db/charon_small.tab` (64 files, 0.89 G minimisers at w41/k19) | **yes**, ONT and Illumina | parameter tuning; laptop users |
| **large** | full reference (old build: 196 files, 16.5 G minimisers) | no: **index ≲ 64 GB** on a server | evaluations for publication |

Tune on the small index, then confirm the chosen settings on the large one. The
laptop numbers below (§5.3, K3) apply to the **small** release. The large release
has a target of about 64 GB for the index. With balanced bins: ONT w41 at FPR 0.01 ≈ 25.5 GB and Illumina
w23 at FPR 0.05 ≈ 54 GB both fit. The current layout (68 GB for ONT w41) doesn't, so balancing is
needed for both releases.

Only two *categories* are needed (host, microbial), but we want to keep many *bins*
because per-bin scoring is believed to improve sensitivity (§3).

## 2. Where we are

Numbers from the March 2025 full build (`dev/index_host_micro4_060325`, w=41, k=19,
FPR 0.01, 3 hashes, 196 files):

| | bins | distinct minimisers | largest bin |
|---|---|---|---|
| human (HPRC, pooled per chromosome) | 25 | 0.26 G (1.6 %) | 21 M |
| bacteria | 153 | 14.3 G | 201 M |
| fungi / archaea / SAR / viruses | 18 | 2.0 G | 224 M |
| **total** | 196 | **16.5 G** | 224 M |

- IBF = 196 bins × 2.77 Gbit = **68 GB uncompressed**. Stored and *queried* in
  seqan3's compressed layout (~39 GB gzipped on disk).
- Only **38 %** of the IBF's capacity is used (16.5 G stored / (196 × 224 M) slots), because
  every bin must be as tall as the largest one (§4.1).
- Build: ~11 h extracting minimisers (20 threads), ~3.5 h filling the IBF (effectively
  single-threaded: every insert is inside one `omp critical`).
- w=31, k=21 doubles the minimiser count (32.4 G), consistent with density
  ≈ 2/(w−k+2). Both of those runs died after creating an ~85 GB IBF.
- Logs roll over at 1 MB × 5 files, so long builds lose their early lines.

## 3. How classification uses bins (and why many bins may matter)

For each read minimiser the IBF returns a bit per bin. Per read, `ReadEntry::get_counts`:

1. totals hits per bin;
2. for each category, picks the **single best bin** (most hits);
3. counts a minimiser as *unique* to a category if it hits that category's best bin and
   not the other category's best bin.

The model is then fitted on unique proportions.

The best-bin step behaves like "compare against the nearest reference group". A
human read is compared only with the *one* microbial bin it resembles most. Minimisers
it happens to share with *other* microbial bins (low-complexity sequence, false
positives spread over 171 bins) don't reduce its human-unique count. If microbes were
collapsed into one filter, every such shared minimiser would count against the read.
That is a real argument for keeping many bins, and it is testable (§7, E2). It changes
my earlier suggestion of one filter per category: that would only be safe if E2 shows
no loss.

## 4. Cost model

### 4.1 Size

```
IBF bits = (number of bins) × (bin height)
bin height = bits_per_element(FPR, h) × (distinct minimisers in the LARGEST bin)
```

| FPR | h=2 | h=3 |
|---|---|---|
| 0.01 (now) | 19.0 | **12.4** |
| 0.02 | 13.1 | 9.5 |
| 0.05 | 7.9 | 6.5 |
| 0.10 | 5.3 | 4.8 |

*(bits per element)*

The lower bound, with perfectly balanced bins, is
`total distinct minimisers × bits_per_element`. For the 196-file build at FPR 0.01
that is 25.5 GB. The current 68 GB is 2.6× that bound because the bins are unbalanced.

### 4.2 Lookup speed: why one IBF beats many separate filters

An IBF query for one minimiser computes h hash positions. At each position it reads one
contiguous **row** holding one bit per bin. With ≤ 512 bins a row fits in a single
64-byte cache line, so a query costs about **h cache misses whether there are 2 bins or
500**.

The same bins as B separate Bloom filters cost about h × B cache misses. So the concern
is right: separate filters only make sense for a *handful* of components (e.g. one
IBF per category = 2h misses). Many bins should stay interleaved.

The compressed layout (`sdsl::sd_vector`) turns each row read into a select/rank
operation, which is much slower than a plain memory read. An uncompressed IBF that
fits in RAM is the fast option.

### 4.3 Per-read work outside the IBF (doesn't change results)

- `get_compression_ratio` gzips every read ([utils.cpp](../../src/utils.cpp)).
- `ReadEntry::bits_` stores a heap-allocated bitvector **per minimiser**.
- `max_bits_` is built for every read but only used by the uncalled debug printer
  `print_result`. It is also wrong: `resize` followed by `emplace_back` doubles its length.
- One `omp critical` per read when adding results.

## 5. Options

Each option lists its effect on size (S), speed (T), accuracy risk (A) and mergeability (M).

### 5.1 Bin layout (removing imbalance) — keeps many bins

**L1. Balance the inputs.** The new `db/` already chunks microbes into ~100 MB
gzipped files. But compressed bytes only roughly predict *distinct minimisers*: redundant
strains (e.g. `bacteria_11_Escherichia_coli`) compress well *and* have few distinct
minimisers. Human chr1 vs chr21 is still about 5×.
S: partial. T: none. A: none. Effort: low (measure, then re-chunk by minimiser count).

**L2. Split large bins into fixed-height technical bins by hash partition.**
Choose a bin height H. A reference group whose minimisers exceed H's capacity is
split into s technical bins, with minimiser x going to partition `x mod s`.
At query time, for that group only the bit of partition `x mod s` is read. It is in the
**same row/cache line**, so there is no extra lookup and no extra false-positive rate.
The group's hit is exactly what one tall bin would give. Small groups can share a
technical bin only if we accept losing them as separate groups (same as `--optimize`).
S: removes most of the imbalance. T: none (slightly more bins per row). A: none in
principle (verify). M: enables §5.6. Effort: medium (layout + bit remapping in
`get_counts`; the bin type grows from `uint8_t`).

**L3. One IBF per category, each with its own height.** Human bins (~10 M) stop paying
for the fungi bin's 224 M height. S: large for human-heavy waste, none within microbes.
T: 2h cache misses instead of h. A: none. Effort: low–medium. Combines with L2.

**L4. Hierarchical IBF (HIBF).** The published general solution to unbalanced bins:
large user bins split, small ones merged, arranged in a hierarchy
(Mehringer et al. 2023, *Genome Biology*, doi:10.1186/s13059-023-02971-4; library
`seqan/hibf`). It scales to thousands of bins. S/T: good. A: should be none. Effort:
high (new dependency, different query API). It's worth knowing about, but L2 + L3 may
get most of the benefit with far less change.

### 5.2 Filter parameters

**P1. FPR / number of hashes** (table in §4.1). Not currently exposed on the command line
(`max_fpr`, `num_hash`). Going from FPR 0.01 to 0.05 is about 2× smaller. h=2 instead of
3 saves one cache miss per query but costs bits. A: medium. The model fits per run, so
a higher background rate may be absorbed. Measure.

### 5.3 Minimiser scheme

**K1. Window w and k.** Density ≈ 2/(w−k+2).

| w, k | minimisers per 500 bp ONT read | per 150 bp Illumina read | index size vs now |
|---|---|---|---|
| 19, 19 (all k-mers) | 482 | 132 | 12× |
| 23, 19 | 161 | 44 | 4× |
| 31, 19 | 69 | 19 | 1.7× |
| 41, 19 (now) | 40 | 11 | 1× |
| 51, 19 | 28 | 8 | 0.71× |
| 61, 19 | 22 | 6 | 0.55× |

**Illumina needs a small w** (11 minimisers per read is why results were poor). Since
each technology can have its own index, the ONT index keeps w ≈ 41 (or larger only if
needed), and the Illumina index gets a smaller w.

The catch is size. With balanced bins (lower bound, h = 3), the full 196-file reference
needs:

| w, k | size vs w41 | FPR 0.01 | FPR 0.05 | FPR 0.10 |
|---|---|---|---|---|
| 41, 19 (ONT now) | 1× | 25.5 GB | 13.5 GB | 9.9 GB |
| 31, 19 | 1.7× | 43.8 | 23.1 | 17.0 |
| 27, 19 | 2.4× | 61.3 | 32.4 | 23.8 |
| 23, 19 | 4× | 102 | 54.0 | 39.7 |

So on the full reference **the Illumina index is the harder one to fit in 16 GB**.
Getting there depends on reducing the reference (R1, dereplication) much more than on
filter tuning. k can also differ: Illumina reads are accurate enough for a longer k
(fewer chance matches), although that doesn't reduce the number of minimisers on its own.

**K2. Query with a smaller w than the index — largely ruled out.** Querying with
w_q < w_i still only finds minimisers that *are in the index*, and a read contains about
`L × density(w_i)` of those whatever w_q is. The only gain is a few extra index minimisers
whose windows run past the read ends. So K2 can't give a w=41 index more signal for
150 bp reads.

**K3. Illumina within 16 GB (small release).** Small index, 0.89 G minimisers at
w41/k19; other (w, k) scaled by density, which slightly overestimates at small w
(repeated k-mers are counted once):

| w, k | minimisers | balanced, FPR 0.01 | 0.05 | 0.10 | current layout (43 % full), 0.01 | per 2×150 bp pair |
|---|---|---|---|---|---|---|
| 41, 19 | 0.89 G | 1.4 GB | 0.7 | 0.5 | 3.2 GB (measured) | 22 |
| 31, 19 | 1.5 G | 2.4 | 1.2 | 0.9 | 5.5 | 37 |
| 27, 19 | 2.1 G | 3.3 | 1.7 | 1.3 | 7.7 | 52 |
| 23, 19 | 3.6 G | 5.5 | 2.9 | 2.1 | 12.8 | 88 |
| 25, 21 | 3.6 G | 5.5 | 2.9 | 2.1 | 12.8 | 86 |
| 31, 25 | 2.7 G | 4.1 | 2.2 | 1.6 | 9.6 | 63 |
| 21, 19 | 5.3 G | 8.2 | 4.3 | 3.2 | 19.2 | 132 |

So for the **small** release a small-w Illumina index fits comfortably once bins are
balanced: w = 23 at FPR 0.05 is about 3 GB. Without balancing, the current layout already
puts w = 23 at 12.8 GB, which is the main reason layout work (L2/L3) comes before the w sweep.

For the **large** release (16.5 G minimisers at w41) the same settings are about 18× bigger
(w = 23, FPR 0.05 ≈ 54 GB balanced), so it needs a server. There, too, short-read accuracy
leans on the scoring side. Levers:

- **Thresholds for short reads.** The defaults (`--confidence 7` unique-hit difference,
  `--min_proportion_diff 0.04`) were tuned for long reads. With about 22 minimisers per pair, a
  7-hit margin rejects most reads. The old poor Illumina results may be at least partly
  due to thresholds rather than w. This is cheap to test: no rebuild, just re-score.
- **Use both mates fully.** Already done (both mates' minimisers are pooled per pair).
- **Higher FPR is less costly for Illumina** (§8.2), but at FPR 0.2 one in five absent
  minimisers looks present, which is significant noise against about 13 true minimisers per read.
  It needs measuring.
- **More compact filters (option F1).** Static xor / binary-fuse style filters need about
  1.1–1.25 × log2(1/p) bits per element (about 5 bits at p = 0.05, versus 6.5 for a 3-hash
  Bloom filter). Gain about 20–30 %, but they're static (rebuild to add), which conflicts with
  §5.6. Keep as a later option.
- **Larger k** reduces chance matches (the noise floor), which helps more at high FPR. It
  doesn't change the density for a given w − k.

### 5.4 Reference content

**R1. Dereplicate microbial references** — *ruled out (§8.3); kept for the record.* One representative per species or ANI
cluster before chunking. Bacteria are 87 % of minimisers, and many genera are heavily
redundant. S: potentially large. T: none. A: low for dehosting (needs measuring),
higher if per-strain resolution were ever needed. Effort: data pipeline, not Charon code.
First measure distinct minimisers per chunk and the overlap between chunks.

**R2. Leave the human pangenome as is.** Pooling HPRC haplotypes per chromosome costs
almost nothing (≈ one genome of minimisers). The small index uses CHM13 + GRCh38.

### 5.5 Representation

**Q1. Uncompressed IBF in RAM** once it fits in ~8 GB. Fastest query. The compressed
layout can stay as the on-disk format, decompressed at load.
**Q2. Memory-map an uncompressed index from SSD.** Avoids load time, but random
access over an index larger than free RAM will thrash. Only useful if it nearly fits.

### 5.6 Mergeable / extendable indexes

The IBF's equal-height constraint, which causes the imbalance, is also what makes
merging possible. If **bin height H is fixed** (by FPR and a chosen capacity), then:

- **append bins** = add columns. seqan3 supports `increase_bin_number_to` on
  uncompressed IBFs (a linear re-interleave, with no re-hashing);
- **merge two indexes** with the same (k, w, h, H) = concatenate columns + merge
  metadata;
- **add sequences to a bin** = insert. Store a HyperLogLog sketch per bin (a few KB,
  mergeable) so Charon knows each bin's true count and can warn when the FPR target is exceeded;
- **large references** fit by being split across several technical bins (L2);
- **removal** isn't possible in a Bloom filter. Rebuild the affected bins instead.

With L3, each category is its own fixed-height IBF, so host and microbial parts can be
built, rebuilt or swapped (e.g. a different host) independently.

### 5.7 Building on the laptop

- The minimiser set per file (up to ~430 M × 8 B plus hash-table overhead) **per thread**
  doesn't fit 8 threads in 16 GB. Bloom insertion is idempotent, so de-duplication isn't
  needed for the filter itself. Stream minimisers straight into the IBF and size bins
  with a HyperLogLog pass (or from bases × density).
- Parallel fill: give each thread a disjoint **range of rows** (positions). Every
  thread hashes, but only writes inside its own range, so no lock is needed.
  Alternatively, insert with atomic OR on 64-bit words.
- No temporary `.min` files. That's ~130 GB of disk traffic for the full build.

## 6. Evaluation setup

### 6.1 Reads

Hall et al. "artificial real" metagenomes: real reads with known source,
Zenodo 10.5281/zenodo.10472796. From Hall MB *et al.*, "Pangenome databases improve
host removal and mycobacteria classification from clinical metagenomic data",
*GigaScience* 2024, doi:10.1093/gigascience/giae010. Stored in
`Projects/Charon/evaluation/hall2024_zenodo_10472796/`, with md5 checksums verified.

| | ONT (`all.ont.fq.gz`) | Illumina (`all_{1,2}.illumina.fq.gz`) |
|---|---|---|
| human | 141,951 reads (3 × 1 Gbp, 1000G samples, R10.4) | NovaSeq, same samples |
| M. tuberculosis | 101,500 reads (R10.3) | HiSeq 4000 |
| Zymo D6322 (7 bacteria + S. cerevisiae) | 813,709 reads (R10.4) | MiSeq |

For development, use fixed stratified subsamples (e.g. 50 k ONT reads, 500 k Illumina
pairs) so one run takes minutes. Confirm the chosen settings on the full sets.

The sister dataset (Zenodo 8339789, simulated, includes viruses) can be added later.

### 6.2 Dev index

`db/charon_small.tab`: 64 bins.
- human: CHM13 + GRCh38, one bin per chromosome (25). The Hall human samples are
  **not** in the reference, which is a fair test.
- microbial: 39 chunks of ~100 MB (17 bacteria, 11 fungi, 2 archaea, 5 SAR, 4 viruses).
  M. tuberculosis and all eight Zymo species **are** present (`bacteria_01, 02, 03, 04,
  10–12, 14, 17`, `fungi_05`), so microbial recall here is optimistic. A
  leave-one-chunk-out run would test reads whose genomes are missing from the reference.

### 6.3 Metrics (per dataset, per read-length stratum)

- **Host leakage**: fraction of human reads *not* called host (they'd end up in a released
  microbial file). This is the safety metric.
- **Microbial loss**: fraction of microbial reads called host.
- Unclassified rate per class.
- Stratify ONT by read length (500–1 k, 1–5 k, 5–20 k, > 20 k bp).
- **Resources**: index size (disk, in RAM), build wall time and peak RSS, classify wall
  time and peak RSS (`/usr/bin/time -l`), index load time.

### 6.4 Tooling needed

- `charon index` options: `--fpr`, `--num-hash` (exist internally, not exposed).
- `charon index --stats-only`: minimisers per bin (+ HLL) without building the filter.
  Lets most of the size modelling run in minutes.
- A small script: assignments + `read2truth` → metrics table, one row per run, to compare
  settings side by side.

## 7. Experiments, in order

| # | change | needs rebuild? | question |
|---|---|---|---|
| E0 | baseline: current code, w41 k19 FPR .01, small index | yes | reference numbers, ONT + Illumina |
| E1 | read-side fixes (§4.3) | no | same calls (diff = 0), how much faster? |
| E2 | scoring: best-bin vs per-category OR (compute both from the same bits) | no | does collapsing bins hurt? (§3) |
| E3 | Illumina scoring/threshold sweep on the E0 index (confidence, proportion diff, min hits) | no | how much of the poor Illumina result is thresholds, not w? |
| E4 | `--stats-only` across (w, k) for both technologies | no filter | minimisers per bin, imbalance, overlap between chunks → size model |
| E5 | FPR sweep 0.01 / 0.02 / 0.05, h 2 / 3 | yes | accuracy vs size, per technology |
| E6 | layout L2 (+ L3) and uncompressed IBF | yes | size and speed at fixed accuracy |
| E7 | w/k sweep, run separately for ONT (around w41) and Illumina (w ≈ 19–31) | yes | best (w, k) per technology within budget |
| E8 | streaming / parallel build | yes | build fits 16 GB, 8 threads |

E1–E3 reuse one index and give the most information for the least compute. From E5 on,
every setting is evaluated on **both** datasets, with ONT reads against the ONT index and
Illumina reads against the Illumina index.

## 8. Decisions and open questions

Decided (2026-10-20):

1. **Per-bin resolution is only needed for scoring.** L2 may merge small groups within
   a category.
2. **Separate indexes per technology are fine.** For Illumina, k may be increased
   (fewer chance matches; slightly fewer distinct minimisers in repetitive sequence).
   A **higher FPR is acceptable for Illumina** than for ONT, because there are fewer
   sequencing errors adding noise.
3. **No further dereplication (R1).** The small index already caps representatives per
   species/genus. The aim is for *every species* to be represented, so that new reads
   are never far from a reference. Size reductions must therefore come from layout, FPR
   and (w, k), not from removing species.
4. **The bin limit can be raised** (`uint8_t` → `uint16_t`). The original cap was a
   proxy for RAM and speed. With fixed-height bins, RAM is set by total content rather
   than bin count, and lookups cost one cache line per 512 bins.

5. **Two releases** (2026-10-20). The *small* release's ONT and Illumina indexes must fit
   a 16 GB laptop. The *large* release runs on a server and is used for publication
   evaluations. Tune on small, confirm on large. For small, a small-w Illumina index is
   feasible once bins are balanced (§5.3, K3).

6. **Large release index target ≈ 64 GB** (2026-10-20).

## 9. Results so far (small index, 2026-10-20)

Setup: `db/charon_small.tab` (64 bins), w41/k19, FPR 0.01, h 3. Built on the laptop in
12.6 min (4.5 min extracting minimisers, 8 min single-threaded fill), 5.8 GB peak. The IBF
is 64 × 400 Mbit = 3.2 GB uncompressed, 2.0 GB on disk (compressed), and 43 % full
(0.89 G minimisers; largest bin fungi_03/06 at 32 M; smallest chrM 0.1 M).
Reads: `ont_60k` (20 k each human / mtb / zymo, shuffled) and `ilmn_300k` (100 k pairs each).
Runs without `--extract`, so the model uses its **default** KDE (see bug B3).
Scripts: `scripts/eval/` (subsample, run + score, compare outputs, offline threshold sweep).

**E0 baseline**

| | human → host | human → microbial | mtb → microbial | zymo → microbial | wall (8 threads) | peak RAM |
|---|---|---|---|---|---|---|
| ONT | 88.96 % (10.8 % filtered by `--min_quality 15`) | 0 | 96.78 % | 99.61 % | 184 s | 2.05 GB |
| Illumina | 95.14 % | 0.001 % | 0 % (all filtered, B5) | 97.73 % | 1,696 s | 1.84 GB |

For ONT reads that pass the quality filter, host recall is 99.4–99.9 % in every length stratum.

**E1 read-side performance: per-read output identical to E0 for every read (evidence and calls).**

| change | ONT | Illumina |
|---|---|---|
| E0 | 184 s | 1,696 s |
| E1a: flat per-hash bit storage, drop `max_bits_`, classify outside the `add_read_to_results` lock once the model is ready, no flush per line | 121 s | 440 s |
| E1b: also format output lines outside the `print_result` lock | 98 s | 221 s (repeat: 312 s, 97 s) |

The profile of E0 Illumina showed the main thread waiting on a lock 95 % of the time. Model
evaluation (KDE `exp`) and printing with `std::endl` were done inside one global critical
section. Illumina timings vary 3× between near-identical runs, so other load on the laptop
(e.g. Dropbox syncing) is interfering. Speed comparisons need repeats on an idle machine.

**E2 best-bin vs category-OR scoring (same index):**

| | best bin | category_or |
|---|---|---|
| ONT human → microbial | 0 / 20,000 | 5 / 20,000 |
| ONT microbial recall | 98.20 % | 99.42 % |
| Illumina human → host | 95.14 % | 93.13 % |
| Illumina human → microbial | 1 / 100,000 | 6 / 100,000 |

Collapsing bins raises microbial sensitivity but leaks host reads and lowers host recall,
which confirms §3. **Keep per-bin scoring.** Size savings must come from balancing (L2/L3), not merging.

**E3 Illumina thresholds** (run with `--min_length 50` so the 75 bp mtb mates are scored;
offline sweep reproduces Charon's calls exactly). Confidence is the only threshold that matters
here (`min_proportion_diff` 0–0.04 makes no difference):

| rule | conf | human → host | mtb | zymo | human → microbial |
|---|---|---|---|---|---|
| current paired (`call_category`, B1) | 7 (default) | 95.2 % | 31.3 % | 97.7 % | 0.001 % |
| | 5 | 96.5 % | 69.7 % | 98.7 % | 0.002 % |
| | 3 | 97.7 % | 91.1 % | 99.1 % | 0.013 % |
| intended (`call_host`), host_lo 0.1 | 7 | 95.9 % | 51.0 % | 98.3 % | 0.002 % |
| | 5 | 97.1 % | 82.6 % | 98.9 % | 0.005 % |
| | 3 | 98.1 % | 94.7 % | 99.2 % | 0.024 % |

At w41, Illumina is a recall-vs-leak trade-off set by the confidence threshold: about 10 minimisers per
75 bp-mate pair, about 22 per 150 bp pair. More minimisers per pair (smaller w for the Illumina
index, E7) should move the whole curve. Even at w41, `call_host` with conf 5 beats the
current default on every class except a small rise in leak.

### Bugs found along the way

| id | where | effect | status |
|---|---|---|---|
| B1 | `process_paired_read_batch` → `add_paired_read(read, r1, r2)` omits `dehost` | paired-end `dehost` uses `call_category` (classify rule), not `call_host` | open |
| B2 | `StatsModel` ctors never set `min_hits_` | `call_category`'s min-hits check reads an uninitialised value | open |
| B3 | `Result` reserves `cached_reads_` only with `--extract` | without `--extract` the model is never fitted, and the default KDE is always used | open |
| B4 | `Result::add_read`/`add_paired_read`, cache-full branch | the read that triggers `force_ready()` is never classified or output (1 read per run) | open |
| B5 | `dehost_main`: paired mode sets `min_length = 80` unconditionally | `--min_length` ignored for pairs, and 75 bp mates (Hall mtb) are silently dropped | **fixed** (explicit `--min_length` honoured) |
| B6 | `index_main`: plog rolls over at 1 MB × 5 | long index builds lose early log lines | open |
| B7 | `ReadEntry::get_counts`: `max_bits_` `resize` then `emplace_back` | doubled vectors (output unused) | **fixed** (removed in E1) |
