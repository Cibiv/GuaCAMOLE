
# GuaCAMOLE

**GC-aware species abundance estimation from metagenomic data**

## Overview
GuaCAMOLE estimates and corrects for the GC bias inherent in most metagenomic sequencing libraries. GuaCAMOLE is based on [Bracken](https://ccb.jhu.edu/software/bracken/) and relies on [Kraken2](https://ccb.jhu.edu/software/kraken2/) for read classification; see our [publication](https://www.biorxiv.org/content/10.1101/2024.09.20.614100) for an in-depth description and evaluation of the algorithm. For instructions how to run GuaCAMOLE on your own dataset see below.

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

## Testing

The `demo_data` folder contains a Kraken2 database containing the 19 bacterial species found in the mock community of Tourlousse *et al.* [[1]](#references), already prepared to be used by Bracken and GuaCAMOLE. The folder also contains a 1% subsample of metagenomic sequencing library [SRR12996245](https://www.ncbi.nlm.nih.gov/sra/?term=SRR12996245) representing that mock community, and the Kraken2 output for that subsample. To run GuaCAMOLE on this data, run

```
./SRR12996245.1pct.sh
```

in that folder. The `SRR12996245.1pct.sh` unzips the FASTQ files and Kraken2 results, changes into to the subdirectory `out/`, and runs GuaCAMOLE with

```
guacamole \
	--output SRR12996245.1pct.guaca \
	--kraken_report ../SRR12996245.1pct_report.txt \
	--kraken_file ../SRR12996245.1pct.kraken \
	--read_files ../SRR12996245.1pct_1.fastq ../SRR12996245.1pct_2.fastq \
	--kraken_db ../demo_db \
	--read_len 150 \
	--fragment_len 400 \
	--length_correction True \
	--threshold 5 \
	--plot True
```

### Docker

GuaCAMOLE is also available as docker image under `laurenz0908/guacamole:latest`. To use the docker image for testing GuaCAMOLE do:

```
mkdir guacamole_out
```

for creating the output folder. Then run the docker image interactively by running

```
docker run -it \
  -v "$(pwd)/guacamole_output:/app/demo_data/out" \
  -w "/app/demo_data/out" \
  laurenz0908/guacamole:latest
```

Then from the docker image shell, run

```
guacamole \
        --output SRR12996245.1pct.guaca \
        --kraken_report ../SRR12996245.1pct_report.txt \
        --kraken_file ../SRR12996245.1pct.kraken \
        --read_files ../SRR12996245.1pct_1.fastq ../SRR12996245.1pct_2.fastq \
        --kraken_db ../demo_db \
        --read_len 150 \
        --fragment_len 400 \
        --length_correction True \
        --threshold 5 \
        --plot True
```

You can exit the image via `exit`. The output files of GuaCAMOLE should now be in the `guacamole_output` directory.

## Usage

### 1. Building a Kraken2 database

To build a standard Kraken2 database compatible with GuaCAMOLE, run

```bash
kraken2-build --fast-build --standard --db path/to/kraken_db \
              --threads number_of_cpus_to_use --no-masking
```

Existing Kraken2 databases can be used, provided that they have been built with the ```--no-masking``` option. Masked databases cannot be used with GuaCAMOLE.

### 2. Create Reference Distribution

GuaCAMOLE requires a Bracken database and a reference distribution, both of which must be created for a read length matching that of the data. By also specifying the ```--fragment_len``` parameter when building the reference distribution, the GC content is computed on a per-fragment instead of a per-read level. This can help with the accuracy of GuaCAMOLE, especially if the fragments are a lot longer than the reads. These databases must be built once for every read length (and fragment length if specified).

```bash
bracken-build -d path/to/kraken_db -t number_of_cpus_to_use -l read_length_of_data 
create-reference-dist --lib_path path/to/kraken_db --ncores number_of_cpus_to_use \
                      --read_len read_length_of_data --fragment_len fragment_length_of_data
```

### 3. Run GuaCAMOLE for Species Abundance Estimation

To estimate species abundances from your data, the reads are first be assigned to taxa with Kraken2, and the Kraken2 output is then processed with GuaCAMOLE. GuaCAMOLE includes Bracken, so no separated invocation of Bracken is necessary.

```bash
kraken2 --db path/to/kraken_db --threads number_of_cpus_to_use --report path/to/kraken_report \
        --paired path/to/reads_1.fastq path/to/reads_2.fastq \
        > path/to/kraken_file
guacamole --kraken_report path/to/kraken_report --kraken_file path/to/kraken_file --kraken_db path/to/kraken_db \
          --read_files path/to/reads_1.fastq path/to/reads_2.fastq
          --read_len read_length_of_data --fragment_len fragment_length_of_data \
          --output result.txt 
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

## References

[1] Tourlousse, D.M., Narita, K., Miura, T. et al, 2021. Validation and standardization of DNA extraction and library construction methods for metagenomics-based human fecal microbiome measurements. *Microbiome* **9**, 95. [https://doi.org/10.1186/s40168-021-01048-3](DOI)
