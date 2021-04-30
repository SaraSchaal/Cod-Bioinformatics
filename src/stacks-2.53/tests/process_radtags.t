#!/usr/bin/env bash

# Preamble
test_path=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
test_data_path="$test_path/"$(basename "${BASH_SOURCE[0]}" | sed -e 's@\.t$@@')
source $test_path/setup.sh

# Setup
barcodes=$test_data_path/frequent_data/Barcodes.txt
freq_in=$test_data_path/frequent_data/in.fastq.gz 
freq_in2=$test_data_path/frequent_data/in.fastq

plan 15

# # Example libtap tests.  Uncomment to run.
# ok "This test will pass" true
# ok "This test will fail" ls -al /this/file/does/not/exist
# diag 'I just love word plays ...'
# ok "This test is expected to fail# TODO fix this" ls -al /neither/does/this/file
# skip "This command doesn't make sense:" more /dev/null

# process_radtags tests

ok_ 'input gzfastq' \
    000_input_gzfastq \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred33 -e sbfI"

ok_ 'input fastq' \
    001_input_fastq \
    "process_radtags -i fastq -f $freq_in2 -o %out -b $barcodes -E phred33 -e sbfI"

diag 'FIXME: Input files for this test are NOT actaully phred64 encoded! This is just an example test...'

ok_ 'input phred64' \
    002_input_phred64 \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred64 -e sbfI"

ok_ 'clean' \
    003_clean_data \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred33 -e sbfI -c"

ok_ 'discarded reads' \
    004_discarded_reads \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred33 -e sbfI -D"

ok_ 'fasta output' \
    005_output_fasta \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred33 -e sbfI -y fasta"

ok_ 'discard low quality reads' \
    006_discard_lq \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred33 -e sbfI -q"

ok_ 'rescue barcodes and radtags' \
    007_rescue_bcrt \
    "process_radtags -i gzfastq -p %in -o %out -b $barcodes -E phred33 -e sbfI -r"

ok_ 'truncate final read length' \
    008_truncate \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred33 -e sbfI -t 50"

ok_ 'set window length' \
    009_winlen \
    "process_radtags -i gzfastq -p %in -o %out -b $barcodes -E phred33 -e sbfI -q -w .12"

ok_ 'minimum window score' \
    010_minscore \
    "process_radtags -i gzfastq -f $freq_in -o %out -b $barcodes -E phred33 -e sbfI -q -s 15"

ok_ 'merge output' \
    011_merge \
    "process_radtags -i gzfastq -f %in/in.fastq.gz -o %out -E phred33 -e sbfI --merge"

ok_ 'Remove sequences marked by Illumina as failing chastity/purity filter' \
    012_filt_ill \
    "process_radtags -i gzfastq -f %in/in.fastq.gz -o %out -E phred33 -e sbfI -b $barcodes --filter_illumina"

ok_ 'Disable checking for RAD site' \
    013_disrc \
    "process_radtags -i gzfastq -f %in/in.fastq.gz -o %out -E phred33 -e sbfI -b $barcodes --disable_rad_check"

ok_ 'Provide distance between barcodes for rescue' \
    014_bcdist \
    "process_radtags -i gzfastq -f %in/in.fastq.gz -o %out -E phred33 -e sbfI -b $barcodes --barcode_dist 1 -r"

# I'm not sure yet what finish() does.
finish
