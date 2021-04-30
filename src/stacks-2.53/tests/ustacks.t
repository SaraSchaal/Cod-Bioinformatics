#!/usr/bin/env bash

# Preamble
test_path=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
test_data_path="$test_path/"$(basename "${BASH_SOURCE[0]}" | sed -e 's@\.t$@@')
source $test_path/setup.sh
freq_in=$test_data_path/frequent_inputs/in.fastq.gz
freq_in2=$test_data_path/frequent_inputs/in2.fastq.gz
freq_in3=$test_data_path/frequent_inputs/in.fasta
freq_in4=$test_data_path/frequent_inputs/in.fasta.gz
freq_in5=$test_data_path/frequent_inputs/in.fastq

plan 21

ref_map

ok_ 'input gzipped fastq' \
    000_ingzfastq \
    "ustacks -t gzfastq -f $freq_in -o %out" \

ok_ 'input fastq' \
    001_infastq \
    "ustacks -t fastq -f $freq_in5 -o %out"

ok_ 'input fasta' \
    002_infasta  \
    "ustacks -t fasta -f $freq_in3 -o %out"

ok_ 'inpus gzipped fasta' \
    003_ingzfasta \
    "ustacks -t gzfasta -f $freq_in4 -o %out"

ok_ 'set sample ID=1 (AKA MySQL column 2)' \
    004_sqlid \
    "ustacks -t gzfastq -f $freq_in -o %out -i 1"

ok_ 'specify minimum depth of coverage required to call a stack' \
    005_mindepcov \
    "ustacks -t gzfastq -f $freq_in -o %out -m 3"

ok_ 'specify maximum distance between stacks' \
    006_maxdistbtw \
    "ustacks -t gzfastq -f $freq_in -o %out -M 3"

ok_ 'specify max distance to align secondary reads to primary stacks' \
    007_maxdist_srps \
    "ustacks -t gzfastq -f $freq_in -o %out -N 5"

skip_ 'retain unused reads' \
    008_retain \
    "ustacks -t gzfastq -f %in -o %out -R"

ok_ 'disable haplotype calling from secondary reads' \
    009_hapcall \
    "ustacks -t gzfastq -f $freq_in -o %out -H"

ok_ 'remove highly repetative (likely error) reads' \
    010_remrep \
    "ustacks -t gzfastq -f %in/in.fastq.gz -o %out -r"

ok_ 'enable deleveraging algorithm' \
    011_deleverage \
    "ustacks -t gzfastq -f %in/in.fastq.gz -o %out -d" \
    zip_test

ok_ 'specify max number of stacks at a de novo locus' \
    012_maxlocus \
    "ustacks -t gzfastq -f %in/in.fastq.gz -o %out --max_locus_stacks 4"

ok_ 'specify chi square significance level for calling heteroz/homozygote' \
    013_alpha0.1 \
    "ustacks -t gzfastq -f %in/in.fastq.gz -o %out --alpha 0.1"

ok_ 'specify chi square significance level for calling heteroz/homozygote' \
    014_alpha0.05 \
    "ustacks -t gzfastq -f %in/in.fastq.gz -o %out --alpha 0.05"

ok_ 'specify chi square significance level for calling heteroz/homozygote' \
    015_alpha0.01 \
    "ustacks -t gzfastq -f %in/in.fastq.gz -o %out --alpha 0.01"

ok_ 'specify chi square significance level for calling heteroz/homozygote' \
    016_alpha0.001 \
    "ustacks -t gzfastq -f %in/in.fastq.gz -o %out --alpha 0.001"
    
ok_ 'For bounded model, specify upper bound' \
    017_bound_high \
    "ustacks -t gzfastq -f $freq_in2 -o %out --model_type bounded --bound_high 0.01"

ok_ 'For bounded model, specify lower bound' \
    018_bound_low \
    "ustacks -t gzfastq -f $freq_in2 -o %out --model_type bounded --bound_low 0.001"

ok_ 'For bounded model, specify upper and lower bounds' \
    019_bounded \
    "ustacks -t gzfastq -f $freq_in2 -o %out --model_type bounded --bound_low 0.001 --bound_high 0.01"

ok_ 'For fixed model, specify barcode error frequency rate' \
    020_fixed \
    "ustacks -t gzfastq -f %in/in.fastq.gz -o %out --model_type fixed --bc_err_freq 0.98"

finish