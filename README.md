# Charon
_Clean Host Associated Reads Out Nanopore_

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

### Download an Index

Alternatively, a pre-built index is available if you know a free way to make a 4.1GB file publicly available. 

This index includes 10 references from [HPRC](https://humanpangenome.org/) and representatives of Bacteria, Viruses, Archaea and Fungi with complete genomes downloaded from NCBI RefSeq.
The category names are `[microbial, human]`. 

The uncompressed index has size approximately 6GB and needs to be decompressed with bgzip before use.

### Classify

Classify `reads.fq.gz` and extract the microbial fraction:

```
charon classify -t 8 --db <example.tab.idx> <reads.fq.gz>
```

Additionally extract the microbial fraction of the input dataset

```
charon classify -t 8 --db <example.tab.idx> <reads.fq.gz> --extract microbial --extract_file <reads.microbial.fq.gz>
```

## Installation

Currently available to build from source. 
It has been developed on MacOS and Unix. 
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

### Classify
```
Classify read file using index.
Usage: bin/charon classify [OPTIONS] <fastaq>

Positionals:
<fastaq> FILE [required]    Fasta/q file

Options:
-h,--help                   Print this help message and exit
--chunk_size INT            Read file is read in chunks of this size, to be processed in parallel within a chunk. [default: d]
-t,--threads INT            Maximum number of threads to use. [default: ]
--db FILE [required]        Prefix for the index.
-e,--extract STRING         Reads from this category in the index will be extracted to file.
--extract_file FILE         Fasta/q file for output
--log FILE                  File for log
-v                          Verbosity of logging. Repeat for increased verbosity
```
Outputs:
STDOUT a tab separated file with the following fields
1. `[U,C]` unclassified or classified
2. `read_id`
3. `call` which of the index catagories it has been called as ("" if no call)
4. `length` number of hashed kmers in the read
5. `confidence_score` the defined as log10(highest_prob/second_highest_prob). If you want the winning probability to take the call, set this to 0. For cautious use set this between 10-20. 
6. `details` a space separated breakdown. For each category, lists `category_name:count_hits:proportion_of_hits_unique_to_category:assigned_probability`. Here the probability score is 
the probability of seeing this combination of unique_proportions if the read comes from this category.

If `extract_file` and `--extract` specified, will output a file with a subset of input reads. 