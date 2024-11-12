
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

To install the quadratic programming solver:

```bash
pip install qpsolvers['cvxopt']
```

## Usage

### 1. Create Reference Distribution

GuaCAMOLE relies on the library.fna files downloaded by Kraken2 to generate the reference distributions. The Kraken2 database needs to be downloaded with the ```--no-masking``` option with the ```kraken2-build``` command!

To create a GC reference distribution from a Kraken2 database run the following command with the parameters of your sample:

```bash

create-reference-dist --lib_path path/to/kraken_db --read_len 150 --fragment_len 400 --ncores 20
```

Specifying the ```fragment_len``` parameters for paired-end sequencing samples can help with the accuracy of GuaCAMOLE, especially if the fragments are a lot longer than the reads. If specified, GuaCAMOLE will use the mean GC content of the read pair as an approximation to the GC content of the fragment. If not specified just single reads.

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
- `--threshold`: Minimum number of reads found for a species to estimate its abundance, default=500
- `--length_correction`, genome size correction (taxonomic vs. read abundance), default=False
- `--plot`, True if detailed plots should be generated
- `--fp_cycles`, Number of iterations for false positive removal, default=4
- `--reg_weight`, Determines how strong the regularization should be [between 0 and 1]', default=0.01
- `--fragment_len`, length of the fragment if paired end and known', default=None
- `--fasta`, True if reads are in fasta format, false if fastq', default=False

### Output

The Output is the same as the tab-delimited Bracken output file. Three additional columns are added:

- `Bracken_estimate`, the abundance estimate from Bracken (if `length_correction=True` they are genome length corrected the same as the GuaCAMOLE estimates)
- `GuaCAMOLE_estimate`, the abundances estimated using the abundance parameter from the GuaCAMOLE algorithm
- `GuaCMAOLE_est_eff`, the abundances computed using the estimated efficiencies by GuaCAMOLE (this does also include esimates for the taxa that were labelled as false positives by GuaCAMOLE)
- `GC_content`, the GC content of the taxon's genome

### Demo

To try out GuaCAMOLE you can download the fastq files of one of the replicates of the sample `GH` sequenced for a publication by Tourlouss et al. (2021) (for more details see their and our manuscripts). To download the data install `sra-tools` and run:

```bash
fasterq-dump SRR12996167
```

You need to have a kraken2 database installed with the `--no-masking` flag set. Download the Kraken2 database using:

```bash
kraken2-build --standard --db demo --threads 24 --no-masking
```

Now you can classify the reads using

```bash
kraken2 --db standard --threads 24 --report SRR12996167_report.txt --paired SRR12996167_1.fastq SRR12996167_2.fastq > SRR12996167.kraken
```

Then build a Bracken database for the read length of the sequencing data which is 150 bp:

```bash
bracken-build -d standard -t 24 -l 150 
```

Now build a GuaCAMOLE database for the corresponding read and fragment length (should take around an hour):

```bash
crate-reference-dist --lib_path standard --read_len 150 --ncores 20 --fragment_len 300
```

Now run GuaCAMOLE using (should take around 20 minutes):

```bash
guacamole --kraken_report SRR12996167_report.txt --kraken_file SRR12996167.kraken --kraken_db standard --read_len 150 fragment_len 150 --length_correction True --output SRR12996167_guacamole.out --read_files SRR12996167_1.fastq SRR12996167_2.fastq
```
