#!/usr/bin/env python

import optparse
import sys
import os
import gzip
import re

#
# Global configuration variables.
#
path     = ""
aln_path = ""
out_path = ""
batch_id = -1

def parse_command_line():
    global path
    global batch_id

    p = optparse.OptionParser()

    #
    # Add options.
    #
    p.add_option("-p", action="store", dest="path",
                 help="path to Stacks directory.")
    p.add_option("-b", action="store", dest="batch_id",
                 help="Stacks batch ID.")

    #
    # Parse the command line
    #
    (opts, args) = p.parse_args()

    if opts.path != None:
        path = opts.path
    if opts.batch_id != None:
        batch_id = int(opts.batch_id)

    if len(path) == 0 or os.path.exists(path) == False:
        print >> sys.stderr, "You must specify a valid path to Stacks input files."
        p.print_help()
        sys.exit()

    if batch_id < 0:
        pritn >> sys.stderr, "You must specify the batch ID that was supplied to Stacks."
        p.print_help()
        sys.exit()

    if path.endswith("/") == False:
        path += "/"
        
        
def find_stacks_files(path, files):
    try:
        entries = os.listdir(path)

        for entry in entries:
            pos = entry.find(".matches.tsv.gz")
            if (pos == -1):
                pos = entry.find(".matches.tsv")
            if (pos != -1):
                files.append(entry[0:pos])
        print >> sys.stderr, "Found", len(files), "Stacks samples."

    except:
        print >> sys.stderr, "Unable to read files from Stacks directory, '" + path + "'"


def parse_cigar(cigar, components):
    #
    # Parse a cigar string, e.g. 48M1D47M or 43M3D52M3D.
    #
    start = 0
    end   = 0
    dist  = ""
    for c in cigar:
        if c.isalpha():
            dist   = int(cigar[start:end])
            op     = c.upper();
            end   += 1
            start  = end
            components.append((op, dist))
        else:
            end += 1


def adjust_snps(cigar, sample_snps):
    #
    # Adjust SNP positions according to how this sample was aligned to the catalog
    # (which is described by the CIGAR string).
    #
    if len(sample_snps) == 0:
        return

    comp = []
    parse_cigar(cigar, comp)

    index  = 0
    offset = 0
    bp     = 0
    for (op, dist) in comp:
        if op == 'M' or op == 'I':
            bp += dist
            while index < len(sample_snps) and sample_snps[index] < bp:
                sample_snps[index] += offset
                index += 1
        if op == 'D':
            offset += dist


def count_sample_snps(path, file, sample_snps):
    matches = {}
    cigars  = {}
    #
    # Open the matches file and load the matches to the catalog.
    #
    p = path + file + ".matches.tsv.gz"
    if os.path.exists(p):
        gzipped = True;
        fh = gzip.open(p, 'rb')
    else:
        gzipped = False;
        fh = open(path + file + ".matches.tsv", "r")

    for line in fh:
        line = line.strip("\n")

        if len(line) == 0 or line[0] == "#":
            continue

        parts = line.split("\t")

        cat_locus    = int(parts[2])
        sample_locus = int(parts[4])
        if len(parts) == 9:
            cigar = parts[8] if len(parts[8]) > 0 else ""
        else:
            cigar = ""

        if sample_locus not in matches:
            matches[sample_locus] = cat_locus
        else:
            if cat_locus != matches[sample_locus]:
                print >> sys.stderr, "Error: sample locus", sample_locus, "matches more than one catalog locus."

        if len(cigar) > 0:
            if sample_locus not in cigars:
                cigars[sample_locus] = cigar
            else:
                if cigar != cigars[sample_locus]:
                    print >> sys.stderr, "Error: sample locus", sample_locus, "has multiple cigar alignments."

    fh.close()

    #
    # Open the SNPs file and record all the SNP positions found in this sample.
    #
    if gzipped:
        fh = gzip.open(path + file + ".snps.tsv.gz", "rb")
    else:
        fh = open(path + file + ".snps.tsv", "r")

    snps = {}

    for line in fh:
        if line[0] == "#":
            continue

        if len(line) == 0:
            continue

        parts        = line.split("\t")
        sample_locus = int(parts[2])
        col          = int(parts[3])
        model        = parts[4]

        if model != "E":
            continue

        if sample_locus not in matches:
            continue

        if sample_locus not in snps:
            snps[sample_locus] = []
        snps[sample_locus].append(col)

    fh.close()

    #
    # Adjust SNP positions according to the gapped alignments recorded in the CIGAR string.
    #
    for sample_locus in snps:
        if sample_locus in cigars:
            adjust_snps(cigars[sample_locus], snps[sample_locus])        

    snp_cnt = 0

    #
    # Transfer this sample's SNPs to the catalog level dictionary.
    #
    for sample_locus in snps:
        cat_locus = matches[sample_locus]
        
        for col in snps[sample_locus]:
            if cat_locus in sample_snps:
                if col in sample_snps[cat_locus]:
                    sample_snps[cat_locus][col] += 1
                    # print >> sys.stderr, "CatLocus:", cat_locus, "; col:", col
                else:
                    sample_snps[cat_locus][col]  = 1
                    snp_cnt  += 1
            else:
                sample_snps[cat_locus]      = {}
                sample_snps[cat_locus][col] = 1
                snp_cnt  += 1
        
    return snp_cnt


def count_catalog_snps(path, batch_id, catalog_snps):
    #
    # Open the tags file and rewrite it with the alignment coordinates.
    #
    p = path + "batch_" + str(batch_id) + ".catalog.snps.tsv.gz"
    if os.path.exists(p):
        gzipped = True;
        fh = gzip.open(path + "batch_" + str(batch_id) + ".catalog.snps.tsv.gz", "rb")
    else:
        gzipped = False;
        fh = open(path + "batch_" + str(batch_id) + ".catalog.snps.tsv", "r")

    snp_cnt = 0

    for line in fh:
        if line[0] == "#":
            continue

        if len(line) == 0:
            continue

        parts     = line.split("\t")
        cat_locus = int(parts[2])
        col       = int(parts[3])
        model     = parts[4]

        if model != "E":
            continue

        snp_cnt += 1

        if cat_locus not in catalog_snps:
            catalog_snps[cat_locus] = []

        catalog_snps[cat_locus].append(col)

    fh.close()

    return snp_cnt

# #                                                                                             # #
# # ------------------------------------------------------------------------------------------- # #
# #                                                                                             # #

parse_command_line()

files        = []
catalog_snps = {}
sample_snps  = {}

#
# Find all the sample files by looking for the matches files in the Stacks directory.
#
find_stacks_files(path, files)

#
# Count all the SNPs identified in the samples.
#
i = 1
for file in files:
    print >> sys.stderr, "Processing file", str(i), "of", len(files), "['" +  file + "']"
    cnt = count_sample_snps(path, file, sample_snps)
    print >> sys.stderr, "  Found", cnt, "heterozygous SNPs in sample."
    i += 1

total_snps = 0
for locus in sample_snps:
    for col in sample_snps[locus]:
        total_snps += 1
print >> sys.stderr, "Found", total_snps, "variable sites across the population."

#
# Count all the SNPs found in the catalog.
#
print >> sys.stderr, "Processing the catalog"
cnt = count_catalog_snps(path, batch_id, catalog_snps)
print >> sys.stderr, "  Found", cnt, "heterozygous SNPs in the catalog."

#
# Count all the SNPs in the catalog but not in any sample: these are the fixed differences cstacks identified.
#
total_snps = 0
fixed_snps = 0

for locus in catalog_snps:
    if locus not in sample_snps:
        continue
    c = {}
    for col in catalog_snps[locus]:
        c[col] = 1
        total_snps += 1
        if col not in sample_snps[locus]:
            # print >> sys.stderr, "Locus:", locus, "Catalog SNPs:", catalog_snps[locus], "Sample SNPs:", sample_snps[locus]
            fixed_snps += 1
    for col in sample_snps[locus]:
        if col not in c:
            print "Locus:", locus, "col:", col, "Catalog SNPs:", catalog_snps[locus], "Sample SNPs:", sample_snps[locus]

print >> sys.stderr, "Found", total_snps, "SNPs across all samples and in the catalog."
print >> sys.stderr, "Found", fixed_snps, "fixed SNPs only in the catalog."
