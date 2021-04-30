#!/usr/bin/env bash

# Preamble
test_path=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
test_data_path="$test_path/"$(basename "${BASH_SOURCE[0]}" | sed -e 's@\.t$@@')
source $test_path/setup.sh

plan 5

#kmer_filter tests
ok_ "filter out rare kmers" \
    000_rare \
    "kmer_filter -f %in/in.fastq -o %out --rare"

ok_ "filter out overly abundant kmers" \
    001_abundant \
    "kmer_filter -f %in/in.fastq -o %out --abundant" 

ok_ "set max kmer frequency for abundance filtering" \
    002_mkf \
    "kmer_filter -f %in/in.fastq -o %out --abundant --max_k_freq 10"

skip_ "input gzfastq" \
    003_ingzfastq \
    "kmer_filter -i gzfastq -f %in/in.fastq.gz -o out"

ok_ "output fastq" \
    004_outfasta \
    "kmer_filter -y fasta -i $freq_in -o %out"