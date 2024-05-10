# varigraph

## Introduction
A fast and resource-efficient genome graph genotyping tool

## Requirements

Please note the following requirements before building and running the software:

* `Linux` operating system
* cmake version `3.12` or higher
* C++ compiler that supports `C++17` or higher, and the `zlib` library installed (we recommend using GCC version `"7.3.0"` or newer) for building `varigraph`

## Installation

**Install via conda**

```shell
conda create -n varigraph
conda activate varigraph
conda install -c duzezhen varigraph
```

**Building on Linux**

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