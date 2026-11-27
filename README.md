# Charon
_Clean Host Associated Reads Out Nanopore_

<img src="./docs/charon_logo.svg" width="400">

Probabilistically identify the host and microbial reads in a metagenomic dataset.

[TOC]: #

# Table of Contents
- [Introduction](#introduction)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)


## Introduction

Charon allows indexing and classification of long reads into a small number of categories (<256 categories with 
references provided by <256 files) with the additional option of extracting one of these categories.

It was developed to classify metagenomic reads as either microbial or host. Unlike other de-hosting methods, it 
calculates a probability score of a read _belonging to a given category and not the other categories._
Other metagenomic de-hosting methods piggy-back off taxnomic classification methods (which give too specific taxonomic 
information) or mapping methods (which give too specific alignment information), neither of which provided probability
scores which could be used to set an acceptable threshold for e.g. releasing metagenomic data into the public domain.

## Quick Start

### Download an Index

Pre-built indexes are published on Zenodo: one record for
[nanopore reads](https://doi.org/10.5281/zenodo.24559698) and one for Illumina reads. `charon fetch` downloads one,
checks its md5 against Zenodo, and decompresses it:

```
charon fetch --list                        # what is available, with download sizes
charon fetch                               # charon_big.nanopore.idx
charon fetch --host-only                   # charon_big_host_only.nanopore.idx, for dehost
charon fetch --tech illumina --size small -o db/
```

- `--tech` is `nanopore` (the default) or `illumina`. `--size big` (the default) has 200 human assemblies,
  `--size small` has 4. Each has human and microbial categories. `--host-only` gives the big index with just the human category; there is no small host-only index.
- The download size is printed before it starts. Big nanopore is about 6.3 GB to download (7.3 GB on disk); big
  Illumina is 28.7 GB on disk. While decompressing, both the `.gz` and the index are on disk.
- The index is written to `-o DIR`, else `$CHARON_DB` if set, else the current directory.
- An interrupted download resumes when you rerun the same command. Rerunning after it finishes does nothing unless
  Zenodo has a newer version of the index; then `--force` replaces it. Where each index came from is recorded in
  `<index>.source.tsv`.
- `--record ID` fetches a specific Zenodo version instead of the latest, and `--keep-gz` keeps the downloaded `.gz`.
- These indexes need Charon v1.1.7 or later; `charon fetch` checks the index format before downloading.

### Build an Index

This takes a tab separated file as input; the first column specifies the path to a reference file, the second column specified the name of the category which they belong to.

```
$ cat example.tab
/path/to/file1.fq.gz    microbial
/path/to/file2.fq.gz    microbial
/path/to/file3.fq.gz    host
```

The index can then be built with:
```
charon index -t 8 --profile nanopore <example.tab>
```

`--profile` chooses parameters for the read technology the index will classify: `nanopore` (w 41, k 19) or
`illumina` (w 27, k 19, a smaller window so short reads get enough minimisers). Build one index per technology you
need. Unless `--prefix` is given, the index is written to `<example.tab>.w<w>k<k>.idx`, so the two profiles give
`example.tab.w41k19.idx` and `example.tab.w27k19.idx`. `-w` and `-k` can also be set directly.

### Dehost

Classify `reads.fq.gz` using the categories in the index (one of which must be "host" or "human"):

```
charon dehost -t 8 --db <example.tab.w41k19.idx> <reads.fq.gz>
```

Additionally extract the microbial fraction of the input dataset

```
charon dehost -t 8 --db <example.tab.w41k19.idx> <reads.fq.gz> --extract microbial --prefix <prefix>
```

#### Host-only index

An index with only the host category (a single `human` or `host` category in the tab file) can also be used with `dehost`.
There is no microbial category to compare against, so each read is called host or not host: a read is called host when the
model supports it coming from the host distribution, and every other read is left uncalled (`U`). This is the closest
equivalent to other host depletion tools, which also only look for the host, and is intended for evaluating how much the
microbial category adds.

```
charon dehost -t 8 --db <host_only.tab.w41k19.idx> <reads.fq.gz> --extract all --prefix <prefix>
```

Reads which are not called host are extracted as `<prefix>_non_host`. `classify` still needs at least two categories.

## Installation

### Docker image
A docker image is hosted on dockerhub and can be pulled using
```
docker pull rmcolq/charon
```

### Building from source
Charon has been developed on MacOS and Unix. 
Requires a compiler for C++20, cmake >= 3.12 and zlib. libcurl is optional: it is needed for `charon fetch`
(on Debian/Ubuntu `apt install libcurl4-openssl-dev`). Without it, cmake warns and Charon builds without the
ability to download indexes; `-DCHARON_FETCH=OFF` builds without it on purpose.

```
git clone https://github.com/rmcolq/charon.git
cd charon
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Usage 

### Dehost
```
Dehost read file into host and other using index.
Usage: bin/charon dehost [OPTIONS] <fastaq>

Positionals:
  <fastaq> FILE [required]              Fasta/q file
  <fastaq> FILE                         Paired Fasta/q file

Options:
  -h,--help                             Print this help message and exit
  
  --db FILE [required]                  Prefix for the index.
  
  -e,--extract STRING                   Reads from this category in the index will be extracted to file (options host, microbial, all).
  --prefix PATH                         Prefix path for output extracted read files
  
  --chunk_size INT                      Read file is read in chunks of this size, to be processed in parallel within a chunk. [default: 100]
  --lo_hi_threshold FLOAT               Threshold used during model fitting stage to decide if read should be used to train lo or hi distribution. [default: 0.15]
  --num_reads_to_fit INT                Number of reads to use to train each distribution in the model. [default: 5000]
  -d,--dist STRING                      Probability distribution to use for modelling. [default: kde]
  
  --min_length INT                      Minimum read length to classify (per mate for paired reads, default 70 when paired). [default: 140]
  --min_quality INT                     Minimum read quality to classify. [default: 15]
  --min_compression FLOAT               Minimum read gzip compression ratio to classify, a measure of how much information is in the read (default 0.3 when paired). [default: 0]
  --confidence INT                      Minimum difference between the top 2 unique hit counts. [default: 7]
  --host_unique_prop_lo_threshold FLOAT Require non-host reads to have unique host proportion below this threshold for classification (default 0.05 when paired). [default: 0.04]
  --min_proportion_diff FLOAT           Minimum difference between the proportion of (non-unique) kmers found in each category (default 0.28 when paired). [default: 0.07]
  --min_probability_diff FLOAT          Minimum difference between the probability found in each category. [default: 0]
  
  --log FILE                            File for log
  -t,--threads INT                      Maximum number of threads to use. [default: ]
  -v                                    Verbosity of logging. Repeat for increased verbosity
```

The threshold defaults depend on the input. A single read file (nanopore) uses the `[default: …]` values. Paired read
files (Illumina) use the "when paired" values for `--min_length`, `--min_compression`,
`--host_unique_prop_lo_threshold` and `--min_proportion_diff`, unless the option is given. The paired defaults were
chosen so that, in a real-read Illumina mixture of human, MTB and Zymo runs, no human pair left called microbial is
human by BLAST; they cost about 0.4 points of microbial pairs there. `--min_length` applies to each mate, so it must
be below the mate length.

Outputs:

1. `[U,C,F]` unclassified, classified, or filtered
2. `read_id`
3. `call` which of the index categories it has been called as ("" if no call)
4. `length` length of read
5. `num_hashes` number of hashed kmers in the read
6. `mean_quality` mean quality of the read
7. `confidence_score` the difference between the number of hits assigned uniquely to called category vs next highest (capped at 255). 
7. `compression` the gzip compression ratio of the read, a measure of the information/complexity of the read.
8. `details` a space separated breakdown. For each category, lists `category_name:count_hits:proportion_of_hits_for_category:proportion_of_hits_unique_to_category:assigned_probability`. 

The probability score is the relative probability of seeing the number of unique hits against this category if the read is truly from the positive distribution, rather than the negative distribution.

Reads removed before classification, because they are shorter than `--min_length` (for a pair, either read) or have a mean quality below `--min_quality`, are not dropped: each has a row with status `F`, so every read of the input appears in the table. For these rows `call` is empty, `num_hashes` and `confidence_score` are 0, `compression` is `NA`, and `details` holds the reason, `filter:too_short` or `filter:low_quality`. `length` is the sum of both reads for a pair. Filtered reads are not extracted.

If `extract_file` and `--prefix` specified, will output a file with a subset of input reads belonging to the names index category. Specifying `--extract all` will generate a file for both host and microbial reads (excludes unclassified, except for a host-only index where the reads which are not host are written as `non_host`).

---

## Credits

**Design & development:** [Rachel Colquhoun](https://github.com/rmcolq)

**Coding:** [Rachel Colquhoun](https://github.com/rmcolq), with collaboration from ELM Qwen-3.5 for release 1.1.0

---

## Funding

Development of Charon has been supported by:

- **Wellcome Trust** — through the [ARTIC Network](https://artic.network/), a project developing rapid genomic surveillance and outbreak response tools for viral pathogens (Discretionary Award 313694/Z/24/Z and Collaborative Award 206298/Z/17/Z).
- **Gates Foundation** — supporting open-source genomic epidemiology tools for global health.

### The ARTIC Network

The [ARTIC Network](https://artic.network/) is a collaborative project focused on real-time genomic surveillance of viral outbreaks. It develops open laboratory protocols, bioinformatics pipelines, and analysis tools — including Charon — to enable rapid phylogenetic analysis during public health emergencies. ARTIC tools have been deployed in responses to Ebola, SARS-CoV-2, and other viral outbreaks worldwide.
