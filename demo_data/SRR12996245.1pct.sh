#!/bin/sh
set -e

test -e SRR12996245.1pct.kraken || bzip2 -d -k SRR12996245.1pct.kraken.bz2
test -e SRR12996245.1pct_1.fastq || bzip2 -d -k SRR12996245.1pct_1.fastq.bz2
test -e SRR12996245.1pct_2.fastq || bzip2 -d -k SRR12996245.1pct_2.fastq.bz2

mkdir out
cd out
exec guacamole \
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
