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
    global out_path
    global aln_path
    global path
    global batch_id

    p = optparse.OptionParser()

    #
    # Add options.
    #
    p.add_option("-o", action="store", dest="out_path",
                 help="write modified Stacks files to this output path.")
    p.add_option("-a", action="store", dest="aln_path",
                 help="SAM file containing catalog locus alignments.")
    p.add_option("-p", action="store", dest="path",
                 help="path to Stacks directory.")
    p.add_option("-b", action="store", dest="batch_id",
                 help="Stacks batch ID.")

    #
    # Parse the command line
    #
    (opts, args) = p.parse_args()

    if opts.out_path != None:
        out_path = opts.out_path
    if opts.aln_path != None:
        aln_path = opts.aln_path
    if opts.path != None:
        path = opts.path
    if opts.batch_id != None:
        batch_id = int(opts.batch_id)

    if len(out_path) == 0 or os.path.exists(out_path) == False:
        print >> sys.stderr, "You must specify a valid path to write files to."
        p.print_help()
        sys.exit()

    if out_path.endswith("/") == False:
        out_path += "/"

    if len(aln_path) == 0 or os.path.exists(aln_path) == False:
        print >> sys.stderr, "You must specify a valid path to a SAM file containing catalog locus alignments."
        p.print_help()
        sys.exit()

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
        print >> sys.stderr, "Unable to read files from Stacks directory, '", path, "'"


def parse_catalog_alignments(aln_path, alns):
    fh = open(aln_path, "r")

    for line in fh:
	line = line.strip("\n")

	if len(line) == 0 or line[0] == "#" or line[0] == "@":
            continue

	parts = line.split("\t")
        locus = int(parts[0])
        chr   = parts[2]
        bp    = int(parts[3])
        flag  = int(parts[1])

        #
        # Check which strand the read is aligned to.
        #
        if flag & 16 > 0:
            alns[locus] = (chr, bp, "-");
        else:
            alns[locus] = (chr, bp, "+");

    fh.close()
    print >> sys.stderr, "Loaded", len(alns), "catalog locus alignments from '", aln_path, "'."


def convert_sample(path, file, out_path, alns):
    matches = {}
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
        matches[sample_locus] = cat_locus

    fh.close()

    #
    # Open the tags file and rewrite it with the alignment coordinates.
    #
    if gzipped:
        fh = gzip.open(path + file + ".tags.tsv.gz", "rb")
    else:
        fh = open(path + file + ".tags.tsv", "r")

    out_fh = open(out_path + file + ".tags.tsv", "w")

    alns_written = {}

    for line in fh:
        if line[0] == "#":
            out_fh.write(line)
            continue

	if len(line) == 0:
            continue

	parts        = line.split("\t")
        sample_locus = int(parts[2])
        read_type    = parts[6]

        if read_type == "consensus":
            if sample_locus not in matches:
                continue;
            cat_locus = matches[sample_locus]
            if cat_locus not in alns:
                continue;

            (chr, bp, strand) = alns[cat_locus]

            if sample_locus in alns_written:
                alns_written[sample_locus] += 1
            else:
                alns_written[sample_locus]  = 1;
                
            buf = "\t".join(parts[0:3]) + "\t" + chr + "\t" + str(bp) + "\t" + strand + "\t" + "\t".join(parts[6:])
            out_fh.write(buf)

        elif sample_locus in alns_written:
            out_fh.write(line)

    fh.close()
    out_fh.close()

    #
    # Open the SNPs, Alleles, and Matches files and rewrite those that had alignment coordinates.
    #
    if gzipped:
        fh = gzip.open(path + file + ".snps.tsv.gz", "rb")
    else:
        fh = open(path + file + ".snps.tsv", "r")

    out_fh = open(out_path + file + ".snps.tsv", "w")

    for line in fh:
        if line[0] == "#":
            out_fh.write(line)
            continue

	if len(line) == 0:
            continue

	parts        = line.split("\t")
        sample_locus = int(parts[2])

        if sample_locus in alns_written:
            out_fh.write(line)

    fh.close()
    out_fh.close()

    if gzipped:
        fh = gzip.open(path + file + ".alleles.tsv.gz", "rb")
    else:
        fh = open(path + file + ".alleles.tsv", "r")

    out_fh = open(out_path + file + ".alleles.tsv", "w")

    for line in fh:
        if line[0] == "#":
            out_fh.write(line)
            continue

	if len(line) == 0:
            continue

	parts        = line.split("\t")
        sample_locus = int(parts[2])

        if sample_locus in alns_written:
            out_fh.write(line)

    fh.close()
    out_fh.close()

    if gzipped:
        fh = gzip.open(path + file + ".matches.tsv.gz", "rb")
    else:
        fh = open(path + file + ".matches.tsv", "r")

    out_fh = open(out_path + file + ".matches.tsv", "w")

    for line in fh:
        if line[0] == "#":
            out_fh.write(line)
            continue

	if len(line) == 0:
            continue

	parts        = line.split("\t")
        sample_locus = int(parts[2])

        if sample_locus in alns_written:
            out_fh.write(line)

    fh.close()
    out_fh.close()

    #
    # If it exists, open the model file and rewrite it with the alignment coordinates.
    #
    if gzipped:
        if os.path.exists(path + file + ".models.tsv.gz") == False:
            return len(alns_written)
    elif os.path.exists(path + file + ".models.tsv") == False:
        return len(alns_written)

    if gzipped:
        fh = gzip.open(path + file + ".models.tsv.gz", "rb")
    else:
        fh = open(path + file + ".models.tsv", "r")

    out_fh = open(out_path + file + ".models.tsv", "w")

    for line in fh:
        if line[0] == "#":
            out_fh.write(line)
            continue

	if len(line) == 0:
            continue

	parts        = line.split("\t")
        sample_locus = int(parts[2])
        read_type    = parts[6]

        if sample_locus in alns_written:
            if read_type == "consensus":
                cat_locus = matches[sample_locus]
                (chr, bp, strand) = alns[cat_locus]
                buf = "\t".join(parts[0:3]) + "\t" + chr + "\t" + str(bp) + "\t" + strand + "\t" + "\t".join(parts[6:])
                out_fh.write(buf)
            else:
                out_fh.write(line)

    fh.close()
    out_fh.close()

    return len(alns_written)


def convert_catalog(path, batch_id, out_path, alns):
    #
    # Open the tags file and rewrite it with the alignment coordinates.
    #
    p = path + "batch_" + str(batch_id) + ".catalog.tags.tsv.gz"
    if os.path.exists(p):
        gzipped = True;
        fh = gzip.open(path + "batch_" + str(batch_id) + ".catalog.tags.tsv.gz", "rb")
    else:
        gzipped = False;
        fh = open(path + "batch_" + str(batch_id) + ".catalog.tags.tsv", "r")

    out_fh = open(out_path + "batch_" + str(batch_id) + ".catalog.tags.tsv", "w")

    alns_written = {}

    for line in fh:
        if line[0] == "#":
            out_fh.write(line)
            continue

	if len(line) == 0:
            continue

	parts     = line.split("\t")
        cat_locus = int(parts[2])

        if cat_locus not in alns:
            continue;

        (chr, bp, strand) = alns[cat_locus]

        if cat_locus in alns_written:
            alns_written[cat_locus] += 1
        else:
            alns_written[cat_locus]  = 1;
               
        buf = "\t".join(parts[0:3]) + "\t" + chr + "\t" + str(bp) + "\t" + strand + "\t" + "\t".join(parts[6:])
        out_fh.write(buf)

    fh.close()
    out_fh.close()

    if gzipped:
        fh = gzip.open(path + "batch_" + str(batch_id) + ".catalog.snps.tsv.gz", "rb")
    else:
        fh = open(path + "batch_" + str(batch_id) + ".catalog.snps.tsv", "r")

    out_fh = open(out_path + "batch_" + str(batch_id) + ".catalog.snps.tsv", "w")

    for line in fh:
        if line[0] == "#":
            out_fh.write(line)
            continue

	if len(line) == 0:
            continue

	parts     = line.split("\t")
        cat_locus = int(parts[2])

        if cat_locus in alns_written:
            out_fh.write(line)

    fh.close()
    out_fh.close()

    if gzipped:
        fh = gzip.open(path + "batch_" + str(batch_id) + ".catalog.alleles.tsv.gz", "rb")
    else:
        fh = open(path + "batch_" + str(batch_id) + ".catalog.alleles.tsv", "r")

    out_fh = open(out_path + "batch_" + str(batch_id) + ".catalog.alleles.tsv", "w")

    for line in fh:
        if line[0] == "#":
            out_fh.write(line)
            continue

	if len(line) == 0:
            continue

	parts     = line.split("\t")
        cat_locus = int(parts[2])

        if cat_locus in alns_written:
            out_fh.write(line)

    fh.close()
    out_fh.close()

    return len(alns_written)


parse_command_line()

files = []
alns  = {}

find_stacks_files(path, files)
parse_catalog_alignments(aln_path, alns)

i  = 1
for file in files:
    print >> sys.stderr, "Processing file", str(i), "of", len(files), "['" +  file + "']"
    cnt = convert_sample(path, file, out_path, alns)
    print >> sys.stderr, "  Added alignments for", cnt, "loci."
    i += 1

#
# Now process the catalog.
#
print >> sys.stderr, "Processing the catalog"
cnt = convert_catalog(path, batch_id, out_path, alns)
print >> sys.stderr, "  Added alignments for", cnt, "catalog loci."
