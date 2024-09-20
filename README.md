
# GuaCAMOLE

**GC-aware species abundance estimation from metagenomic data**

## Overview
GuaCAMOLE is a tool for estimating species abundance from metagenomic data, considering GC content. This package implements various scripts to generate reference distributions, correct for GC bias, and estimate species abundances.

## Features
- Create GC reference distributions from Kraken2 databases.
- Estimate species abundances using GC-aware correction.
- Generate detailed plots for data visualization.

## Installation

You can install the package directly from a folder:

```bash
pip install .
```

Or from GitHub:

```bash
pip install git+https://github.com/Cibiv/GuaCAMOLE.git
```

## Requirements

- Python 3.10 or higher
- [numpy==1.26.4](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [matplotlib](https://matplotlib.org/)
- [seaborn](https://seaborn.pydata.org/)
- [biopython](https://biopython.org/)
- [scipy](https://scipy.org/)
- [qpsolvers](https://github.com/stephane-caron/qpsolvers)

To install the quadratic programming solver:

```bash
pip install qpsolvers['cvxopt']
```

## Usage

### 1. Create Reference Distribution

GuaCAMOLE relies on the library.fna files downloaded by Kraken2 to generate the reference distributions. The Kraken2 database needs to be downloaded with the ```--no-masking``` option with the ```kraken2-build``` command!

To create a GC reference distribution from a Kraken2 database, run:

```bash

create-reference-dist --lib_path path/to/kraken_db --read_len 150 --fragment_len 400 --ncores 20
```

For GuaCAMOLE to run it also requires a Bracken database with the correct read length to exist in the Kraken database folder.

### 2. Run GuaCAMOLE for Species Abundance Estimation

To estimate species abundances from your data, run:

```bash
guacamole --kraken_report path/to/report --kraken_file path/to/kraken --kraken_db path/to/kraken_db --read_len 150 --output result.txt --read_files path/to/reads_1.fastq path/to/reads_2.fastq
```

### Command-line Options

- `--kraken_report`: Kraken2 report file (required)
- `--kraken_file`: Kraken2 file with classifications (required)
- `--kraken_db`: Path to the Kraken2 database (required)
- `--read_len`: Read length (required)
- `--output`: Output file name (required)
- `--read_files`: Path to input read files (required)