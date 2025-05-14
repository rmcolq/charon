# Charon
_Clean Host Associated Reads Out Nanopore_

<img src="./docs/charon_logo.svg" width="400">

Probabilistically identify the host and microbial reads in a metagenomic dataset.

UNDER ACTIVE DEVELOPMENT - this tool may not work out of the box but feel free to try it/make suggestions and watch the repo to be informed of releases!

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

Alternatively, a pre-built index for long reads is available via [Zenodo](https://zenodo.org/records/15398095). 

This index includes references from [HPRC](https://humanpangenome.org/) and representatives of Bacteria, Viruses, Archaea, Fungi and Sar with (mostly complete) genomes downloaded from NCBI RefSeq, including FDA-ARGOS.
The category names are `[microbial, human]`. 

The compressed index has size approximately 39GB and needs to be decompressed with gzip before use.

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
charon index -t 8 <example.tab>
```

### Dehost

Classify `reads.fq.gz` using the categories in the index (one of which must be "host" or "human"):

```
charon dehost -t 8 --db <example.tab.idx> <reads.fq.gz>
```

Additionally extract the microbial fraction of the input dataset

```
charon dehost -t 8 --db <example.tab.idx> <reads.fq.gz> --extract microbial --prefix <prefix>
```

## Installation

### Docker image
A docker image is hosted on dockerhub and can be pulled using
```
docker pull rmcolq/charon
```

### Building from source
Charon has been developed on MacOS and Unix. 
Requires a compiler for C++14 and cmake > 3.9.

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
  
  --min_length INT                      Minimum read length to classify. [default: 140]
  --min_quality INT                     Minimum read quality to classify. [default: 15]
  --min_compression FLOAT               Minimum read gzip compression ratio to classify (a measure of how much information is in the read. [default: 0]
  --confidence INT                      Minimum difference between the top 2 unique hit counts. [default: 2]
  --host_unique_prop_lo_threshold INT   Require non-host reads to have unique host proportion below this threshold for classification. [default: 0.05]
  --min_proportion_diff FLOAT           Minimum difference between the proportion of (non-unique) kmers found in each category. [default: 0.04]
  --min_probability_diff FLOAT          Minimum difference between the probability found in each category. [default: 0]
  
  --log FILE                            File for log
  -t,--threads INT                      Maximum number of threads to use. [default: ]
  -v                                    Verbosity of logging. Repeat for increased verbosity
```
Outputs:

1. `[U,C]` unclassified or classified
2. `read_id`
3. `call` which of the index categories it has been called as ("" if no call)
4. `length` length of read
5. `num_hashes` number of hashed kmers in the read
6. `mean_quality` mean quality of the read
7. `confidence_score` the difference between the number of hits assigned uniquely to called category vs next highest (capped at 255). 
7. `compression` the gzip compression ratio of the read, a measure of the information/complexity of the read.
8. `details` a space separated breakdown. For each category, lists `category_name:count_hits:proportion_of_hits_for_category:proportion_of_hits_unique_to_category:assigned_probability`. 

The probability score is the relative probability of seeing the number of unique hits against this category if the read is truly from the positive distribution, rather than the negative distribution.

If `extract_file` and `--prefix` specified, will output a file with a subset of input reads belonging to the names index category. Specifying `--extract all` will generate a file for both host and microbial reads (excludes unclassified).
