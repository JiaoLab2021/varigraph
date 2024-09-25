# varigraph

[![GitHub last commit](https://img.shields.io/github/last-commit/JiaoLab2021/varigraph.svg?label=Last%20commit&logo=github&style=flat)](https://github.com/JiaoLab2021/varigraph/releases)
[![Build Status](https://github.com/JiaoLab2021/varigraph/actions/workflows/ci.yaml/badge.svg)](https://github.com/JiaoLab2021/varigraph/actions)

## Introduction
A fast and resource-efficient genome graph genotyping tool

## Requirements

Please note the following requirements before building and running the software:

* `Linux` operating system
* cmake version `3.12` or higher
* C++ compiler that supports `C++17` or higher, and the `zlib` library installed (we recommend using GCC version `"7.3.0"` or newer) for building `varigraph`

## Installation

### Install via conda

```shell
conda create -n varigraph
conda activate varigraph
conda install -c duzezhen varigraph
```

### Building on Linux

Use the following script to build the software:

1. First, obtain the source code.

```shell
git clone https://github.com/JiaoLab2021/varigraph.git
cd varigraph
```

2. Next, compile the software and add the current directory to your system's `PATH` environment variable.

```shell
cmake ./
make -j 5
echo 'export PATH="$PATH:'$(pwd)'"' >> ~/.bashrc
source ~/.bashrc
```

## Usage

### Input Files

* Reference Genome
* VCF File of Population Variants
* Sample File:

```shell
# Sample File
sample1 sample1.r1.fq.gz sample1.r2.fq.gz
sample2 sample2.r1.fq.gz sample2.r2.fq.gz
...
sampleN sampleN.r1.fq.gz sampleN.r2.fq.gz
```

Please note that the Sample file must be formatted exactly as shown above, where each sample is listed with its corresponding read files.

### Running

For convenience, let's assume the following file names for the input:

* `refgenome.fa`
* `input.vcf.gz`
* `samples.cfg`

`varigraph` runs in two steps: the first step builds the genome graph, and the second step performs the genotyping. Here is the specific code:

**1. Building the Genome Graph:**

```shell
varigraph construct -r refgenome.fa -v input.vcf.gz --save-graph graph.bin
```

* Adjustment for Tetraploid Genome:
   * If your VCF file involves variants from a tetraploid genome, include the `--vcf-ploidy 4` parameter.

**2. Performing Genotyping:**

```shell
varigraph genotype --load-graph graph.bin -s samples.cfg --use-depth
```

* Adjustments for Genotyping:
   * Homozygous Samples: For homozygous samples, add `-g hom` to improve genotyping accuracy.
   * Tetraploid Samples: If your samples are tetraploid, adjust the `--sample-ploidy 4` parameter.
   * Use `--use-depth` for accurate genotyping regardless of ploidy.

## Note on GPU Acceleration

* GPU Version: varigraph also has a GPU-enabled version for faster computation if your server is equipped with GPUs.

* Usage:
   * Use `--gpu` to specify GPU usage. For example, `--gpu 0` uses GPU 0.
   * Adjust GPU memory usage with `--buffer` parameter. Smaller values consume less GPU memory.