#!/usr/bin/env bash

# Preamble
test_path=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
test_data_path="$test_path/"$(basename "${BASH_SOURCE[0]}" | sed -e 's@\.t$@@')
source $test_path/setup.sh
freq_in=$test_data_path/frequent_input/in.bam
freq_in2=$test_data_path/frequent_input/in2.bam
freq_in3=$test_data_path/frequent_input/in.sam

plan 12

ok_ 'input bam' \
    000_inbam \
    "pstacks -t bam -f $freq_in -o %out"

ok_ 'input sam' \
    001_insam \
    "pstacks -t sam -f $freq_in3 -o %out"

skip_ 'input bowtie - need data input file type' \
    002_inbowtie \
    "pstacks -t bowtie -f %in/??? -o %out"

ok_ 'minimum coverage depth for contig report' \
    003_mincovdepth \
    "pstacks -t bam -f $freq_in -o %out -m 2"

ok_ 'R^2 significance level of 0.1 for calling homozygote/heterozygote' \
    004_alpha0.1 \
    "pstacks -t bam -f %in/in.bam -o %out --alpha 0.1"

ok_ 'R^2 significance level of 0.05 for calling homozygote/heterozygote' \
    005_alpha0.05 \
    "pstacks -t bam -f %in/in.bam -o %out --alpha 0.05"

ok_ 'R^2 significance level of 0.01 for calling homozygote/heterozygote' \
    006_alpha0.01 \
    "pstacks -t bam -f %in/in.bam -o %out --alpha 0.01"

ok_ 'R^2 significance level of 0.001 for calling homozygote/heterozygote' \
    007_alpha0.001 \
    "pstacks -t bam -f %in/in.bam -o %out --alpha 0.001"

ok_ 'For bounded model, specify upper bound' \
    008_bound_high \
    "pstacks -t bam -f $freq_in2 -o %out --model_type bounded --bound_high 0.5"

ok_ 'For bounded model, specify lower bound' \
    009_bound_low \
    "pstacks -t bam -f $freq_in2 -o %out --model_type bounded --bound_low 0.2"

ok_ 'For bounded model, specify upper and lower bounds' \
    010_bounded \
    "pstacks -t bam -f $freq_in2 -o %out --model_type bounded --bound_high 0.5 --bound_low 0.2"

ok_ 'For fixed model, specify barcode error frequency' \
    011_fixed \
    "pstacks -t bam -f %in/in.bam -o %out --model_type fixed --bc_err_freq 0.39"

finish