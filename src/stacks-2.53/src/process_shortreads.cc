// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011-2015, Julian Catchen <jcatchen@illinois.edu>
//
// This file is part of Stacks.
//
// Stacks is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Stacks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
//

//
// process_shortreads -- clean raw reads using a sliding window approach;
//   split reads by barcode if barcodes provided, correct barcodes
//   within one basepair, truncate reads on request.
//

#include "process_shortreads.h"

//
// Global variables to hold command-line options.
//
FileT  in_file_type  = FileT::unknown;
FileT  out_file_type = FileT::unknown;
string in_file;
string in_file_p1;
string in_file_p2;
string in_path_1;
string in_path_2;
string out_path;
string barcode_file;
char  *adapter_1;
char  *adapter_2;
barcodet barcode_type    = null_null;
bool     retain_header   = false;
bool     filter_adapter  = false;
bool     paired          = false;
bool     clean           = false;
bool     quality         = false;
bool     recover         = false;
bool     interleaved     = false;
bool     merge           = false;
bool     discards        = false;
bool     overhang        = true;
bool     matepair        = false;
bool     filter_illumina = false;
bool     trim_reads      = true;
uint     truncate_seq    = 0;
int      barcode_dist_1  = 1;
int      barcode_dist_2  = -1;
double   win_size        = 0.15;
int      score_limit     = 10;
int      len_limit       = 31;
int      num_threads     = 1;

//
// How to shift FASTQ-encoded quality scores from ASCII down to raw scores
//     score = encoded letter - 64; Illumina version 1.3 - 1.5
//     score = encoded letter - 33; Sanger / Illumina version 1.6+
int qual_offset  = 33;

//
// Handle variable-size barcodes.
//
uint min_bc_size_1 = 0;
uint max_bc_size_1 = 0;
uint min_bc_size_2 = 0;
uint max_bc_size_2 = 0;

//
// Kmer data for adapter filtering.
//
int kmer_size = 5;
int distance  = 1;
int adp_1_len = 0;
int adp_2_len = 0;
AdapterHash adp_1_kmers, adp_2_kmers;

int main (int argc, char* argv[]) {
    IF_NDEBUG_TRY

    parse_command_line(argc, argv);

    //
    // If input files are gzipped, output gziped files, unless the user chooses an output type.
    //
    if (out_file_type == FileT::unknown) {
        if (in_file_type == FileT::gzfastq || in_file_type == FileT::bam)
            out_file_type = FileT::gzfastq;
        else
            out_file_type = FileT::fastq;
    }

    cerr << "Using Phred+" << qual_offset << " encoding for quality scores.\n"
         << "Reads trimmed shorter than " << len_limit << " nucleotides will be discarded.\n";
    if (truncate_seq > 0)
        cerr << "Reads will be truncated to " << truncate_seq << "bp\n";
    if (filter_illumina)
        cerr << "Discarding reads marked as 'failed' by Illumina's chastity/purity filters.\n";
    if (filter_adapter) {
        cerr << "Filtering reads for adapter sequence:\n";
        if (adapter_1 != NULL) {
            cerr << "  " << adapter_1 << "\n";
            init_adapter_seq(kmer_size, adapter_1, adp_1_len, adp_1_kmers);
        }
        if (adapter_2 != NULL) {
            cerr << "  " << adapter_2 << "\n";
            init_adapter_seq(kmer_size, adapter_2, adp_2_len, adp_2_kmers);
        }
        cerr << "    " << distance << " mismatches allowed to adapter sequence.\n";
    }

    vector<pair<string, string> >        files;
    vector<BarcodePair>                  barcodes;
    set<string>                          se_bc, pe_bc;
    map<BarcodePair, ofstream *>         pair_1_fhs, pair_2_fhs, rem_1_fhs, rem_2_fhs;
    map<BarcodePair, gzFile *>           pair_1_gzfhs, pair_2_gzfhs, rem_1_gzfhs, rem_2_gzfhs;
    map<string, map<string, long> >      counters;
    map<BarcodePair, map<string, long> > barcode_log;

    build_file_list(files);
    load_barcodes(barcode_file, barcodes, se_bc, pe_bc, min_bc_size_1, max_bc_size_1, min_bc_size_2, max_bc_size_2);
    if (recover && barcode_type != null_null) {
        if (barcode_type == index_null || barcode_type == inline_null)
            cerr << "Will attempt to recover barcodes with at most " << barcode_dist_1 << " mismatches.\n";
        else
            cerr << "Will attempt to recover barcodes with at most " << barcode_dist_1 << " / " << barcode_dist_2 << " mismatches.\n";
    }

    if (out_file_type == FileT::gzfastq || out_file_type == FileT::gzfasta)
        open_files(files, barcodes, pair_1_gzfhs, pair_2_gzfhs, rem_1_gzfhs, rem_2_gzfhs, counters);
    else
        open_files(files, barcodes, pair_1_fhs, pair_2_fhs, rem_1_fhs, rem_2_fhs, counters);

    int result = 1;

    for (uint i = 0; i < files.size(); i++) {
        cerr << "Processing file " << i+1 << " of " << files.size() << " [" << files[i].first.c_str() << "]\n";

        counters[files[i].first]["total"]        = 0;
        counters[files[i].first]["ill_filtered"] = 0;
        counters[files[i].first]["low_quality"]  = 0;
        counters[files[i].first]["trimmed"]      = 0;
        counters[files[i].first]["adapter"]      = 0;
        counters[files[i].first]["ambiguous"]    = 0;
        counters[files[i].first]["retained"]     = 0;
        counters[files[i].first]["orphaned"]     = 0;
        counters[files[i].first]["recovered"]    = 0;

        if (paired) {
            if (out_file_type == FileT::gzfastq || out_file_type == FileT::gzfasta)
                    result = process_paired_reads(files[i].first, files[i].second,
                                              se_bc, pe_bc,
                                              pair_1_gzfhs, pair_2_gzfhs, rem_1_gzfhs, rem_2_gzfhs,
                                              counters[files[i].first], barcode_log);
            else
                result = process_paired_reads(files[i].first, files[i].second,
                                              se_bc, pe_bc,
                                              pair_1_fhs, pair_2_fhs, rem_1_fhs, rem_2_fhs,
                                              counters[files[i].first], barcode_log);
        } else {
            if (out_file_type == FileT::gzfastq || out_file_type == FileT::gzfasta)
                    result = process_reads(files[i].first,
                                       se_bc, pe_bc,
                                       pair_1_gzfhs,
                                       counters[files[i].first], barcode_log);
            else
                result = process_reads(files[i].first,
                                       se_bc, pe_bc,
                                       pair_1_fhs,
                                       counters[files[i].first], barcode_log);
        }
        cerr <<        "  "
             << counters[files[i].first]["total"] << " total reads; ";
        if (filter_illumina)
            cerr << "-" << counters[files[i].first]["ill_filtered"] << " failed Illumina reads; ";
        cerr
             << "-" << counters[files[i].first]["ambiguous"]   << " ambiguous barcodes; "
             << "+" << counters[files[i].first]["recovered"]   << " recovered; "
             << "-" << counters[files[i].first]["low_quality"] << " low quality reads; "
             << counters[files[i].first]["retained"] << " retained reads.\n"
             << "    ";
        if (filter_adapter)
            cerr << counters[files[i].first]["adapter"] << " reads with adapter sequence; ";
        cerr << counters[files[i].first]["trimmed"]     << " trimmed reads; "
             << counters[files[i].first]["orphaned"]    << " orphaned paired-ends.\n";

        if (!result) {
            cerr << "Error processing reads.\n";
            break;
        }
    }

    cerr << "Closing files, flushing buffers...\n";
    if (out_file_type == FileT::gzfastq || out_file_type == FileT::gzfasta) {
        close_file_handles(pair_1_gzfhs);
        if (paired) {
            close_file_handles(rem_1_gzfhs);
            close_file_handles(rem_2_gzfhs);
            close_file_handles(pair_2_gzfhs);
        }
    } else {
        close_file_handles(pair_1_fhs);
        if (paired) {
            close_file_handles(rem_1_fhs);
            close_file_handles(rem_2_fhs);
            close_file_handles(pair_2_fhs);
        }
    }

    print_results(argc, argv, barcodes, counters, barcode_log);

    return 0;
    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

template<typename fhType>
int
process_paired_reads(string prefix_1,
                     string prefix_2,
                     set<string> &se_bc, set<string> &pe_bc,
                     map<BarcodePair, fhType *> &pair_1_fhs,
                     map<BarcodePair, fhType *> &pair_2_fhs,
                     map<BarcodePair, fhType *> &rem_1_fhs,
                     map<BarcodePair, fhType *> &rem_2_fhs,
                     map<string, long> &counter,
                     map<BarcodePair, map<string, long> > &barcode_log) {
    Input *fh_1=NULL, *fh_2=NULL;
    RawRead  *r_1=NULL, *r_2=NULL;
    ofstream *discard_fh_1=NULL, *discard_fh_2=NULL;

    int return_val = 1;

    string path_1 = in_path_1 + prefix_1;
    string path_2 = in_path_2 + prefix_2;

    if (interleaved)
        cerr << "  Reading data from:\n  " << path_1 << "\n";
    else
        cerr << "  Reading data from:\n  " << path_1 << " and\n  " << path_2 << "\n";

    if (in_file_type == FileT::fastq) {
        fh_1 = new Fastq(path_1.c_str());
        fh_2 = interleaved ? fh_1 : new Fastq(path_2.c_str());
    } else if (in_file_type == FileT::gzfastq) {
        fh_1 = new GzFastq(path_1.c_str());
        fh_2 = interleaved ? fh_1 : new GzFastq(path_2.c_str());
    } else if (in_file_type == FileT::bam) {
        fh_1 = new BamUnAln(path_1.c_str());
        fh_2 = fh_1;
    } else if (in_file_type == FileT::bustard) {
        fh_1 = new Bustard(path_1.c_str());
        fh_2 = interleaved ? fh_1 : new Bustard(path_2.c_str());
    }

    //
    // Open a file for recording discarded reads
    //
    if (discards) {
        path_1 = out_path + prefix_1 + ".discards";
        discard_fh_1 = new ofstream(path_1.c_str(), ifstream::out);

        if (discard_fh_1->fail()) {
            cerr << "Error opening discard output file '" << path_1 << "'\n";
            exit(1);
        }

        path_2 = out_path + prefix_2 + ".discards";
        discard_fh_2 = new ofstream(path_2.c_str(), ifstream::out);

        if (discard_fh_1->fail()) {
            cerr << "Error opening discard output file '" << path_2 << "'\n";
            exit(1);
        }
    }

    //
    // Read in the first record, initializing the Seq object s. Then
    // initialize the Read object r, then loop, using the same objects.
    //
    Seq *s_1 = fh_1->next_seq();
    Seq *s_2 = fh_2->next_seq();
    if (s_1 == NULL || s_2 == NULL) {
        cerr << "Attempting to read first pair of input records, unable to allocate "
             << "Seq object (Was the correct input type specified?).\n";
        exit(1);
    }

    r_1 = new RawRead(strlen(s_1->seq), 1, min_bc_size_1, win_size);
    r_2 = new RawRead(strlen(s_2->seq), 2, min_bc_size_2, win_size);

    BarcodePair bc;
    //
    // If no barcodes were specified, set the barcode object to be the input file names.
    //
    if (max_bc_size_1 == 0)
        bc.set(prefix_1, prefix_2);

    long i = 1;

    do {
        if (i % 10000 == 0) cerr << "  Processing short read " << i << "       \r";

        parse_input_record(s_1, r_1);
        parse_input_record(s_2, r_2);
        counter["total"] += 2;

        if (barcode_type != null_null &&
            barcode_type != inline_null &&
            barcode_type != index_null)
            bc.set(r_1->se_bc, r_2->pe_bc);
        else if (barcode_type != null_null)
            bc.set(r_1->se_bc);

        process_barcode(r_1, r_2, bc, pair_1_fhs, se_bc, pe_bc, barcode_log, counter);

        //
        // Adjust the size of the read to accommodate truncating the sequence and variable
        // barcode lengths. With standard Illumina data we want to output constant length
        // reads even as the barcode size may change. Other technologies, like IonTorrent
        // need to be truncated uniformly.
        //
        if (truncate_seq > 0) {
            if (truncate_seq + r_1->inline_bc_len <= r_1->len)
                r_1->set_len(truncate_seq + r_1->inline_bc_len);
            if (truncate_seq + r_2->inline_bc_len <= r_2->len)
                r_2->set_len(truncate_seq + r_2->inline_bc_len);
        } else {
            if (barcode_type == inline_null || barcode_type == inline_inline ||        barcode_type == inline_index)
                r_1->set_len(r_1->len - (max_bc_size_1 - r_1->inline_bc_len));
            if (barcode_type == inline_index ||        barcode_type == index_index)
                r_2->set_len(r_2->len - (max_bc_size_2 - r_2->inline_bc_len));
        }

        if (r_1->retain)
            process_singlet(r_1, false, barcode_log[bc], counter);
        if (r_2->retain)
            process_singlet(r_2, true,  barcode_log[bc], counter);

        if (matepair) {
            rev_complement(r_1->seq, r_1->inline_bc_len, overhang);
            reverse_qual(r_1->phred, r_1->inline_bc_len, overhang);
        }

        int result_1 = 1;
        int result_2 = 1;

        if (r_1->retain && r_2->retain) {
            if (retain_header) {
                result_1 = (out_file_type == FileT::fastq || out_file_type == FileT::gzfastq) ?
                    write_fastq(pair_1_fhs[bc], s_1, r_1) :
                    write_fasta(pair_1_fhs[bc], s_1, r_1);
                result_2 = (out_file_type == FileT::fastq || out_file_type == FileT::gzfastq) ?
                    write_fastq(pair_2_fhs[bc], s_2, r_2) :
                    write_fasta(pair_2_fhs[bc], s_2, r_2);
            } else {
                result_1 = (out_file_type == FileT::fastq || out_file_type == FileT::gzfastq) ?
                    write_fastq(pair_1_fhs[bc], r_1, overhang) :
                    write_fasta(pair_1_fhs[bc], r_1, overhang);
                result_2 = (out_file_type == FileT::fastq || out_file_type == FileT::gzfastq) ?
                    write_fastq(pair_2_fhs[bc], r_2, overhang) :
                    write_fasta(pair_2_fhs[bc], r_2, overhang);
            }

        } else if (r_1->retain && !r_2->retain) {
            //
            // Write to a remainder file.
            //
            if (retain_header)
                result_1 = (out_file_type == FileT::fastq || out_file_type == FileT::gzfastq) ?
                    write_fastq(rem_1_fhs[bc], s_1, r_1) :
                    write_fasta(rem_1_fhs[bc], s_1, r_1);
            else
                result_1 = (out_file_type == FileT::fastq || out_file_type == FileT::gzfastq) ?
                    write_fastq(rem_1_fhs[bc], r_1, overhang) :
                    write_fasta(rem_1_fhs[bc], r_1, overhang);

        } else if (!r_1->retain && r_2->retain) {
            // Write to a remainder file.
            if (retain_header)
                result_2 = (out_file_type == FileT::fastq || out_file_type == FileT::gzfastq) ?
                    write_fastq(rem_2_fhs[bc], s_2, r_2) :
                    write_fasta(rem_2_fhs[bc], s_2, r_2);
            else
                result_2 = (out_file_type == FileT::fastq || out_file_type == FileT::gzfastq) ?
                    write_fastq(rem_2_fhs[bc], r_2, overhang) :
                    write_fasta(rem_2_fhs[bc], r_2, overhang);
        }

        if (!result_1 || !result_2) {
            cerr << "Error writing to output file for '" << bc.str() << "'\n";
            return_val = -1;
            break;
        }

        if (discards && !r_1->retain)
            result_1 = (out_file_type == FileT::fastq || out_file_type == FileT::gzfastq) ?
                write_fastq(discard_fh_1, s_1) :
                write_fasta(discard_fh_1, s_1);
        if (discards && !r_2->retain)
            result_2 = (out_file_type == FileT::fastq || out_file_type == FileT::gzfastq) ?
                write_fastq(discard_fh_2, s_2) :
                write_fasta(discard_fh_2, s_2);

        delete s_1;
        delete s_2;

        if (!result_1 || !result_2) {
            cerr << "Error writing to discard file for '" << bc.str() << "'\n";
            return_val = -1;
            break;
        }

        i++;
    } while ((s_1 = fh_1->next_seq()) != NULL &&
             (s_2 = fh_2->next_seq()) != NULL);

    if (discards) {
        delete discard_fh_1;
        delete discard_fh_2;
    }

    delete fh_1;
    if (interleaved == false) delete fh_2;

    return return_val;
}

template<typename fhType>
int
process_reads(string prefix,
              set<string> &se_bc, set<string> &pe_bc,
              map<BarcodePair, fhType *> &pair_1_fhs,
              map<string, long> &counter,
              map<BarcodePair, map<string, long> > &barcode_log) {
    Input *fh=NULL;
    RawRead  *r;
    ofstream *discard_fh=NULL;

    int return_val = 1;

    string path = in_path_1 + prefix;

    if (in_file_type == FileT::fastq)
        fh = new Fastq(path.c_str());
    else if (in_file_type == FileT::gzfastq)
        fh = new GzFastq(path.c_str());
    else if (in_file_type == FileT::bam)
        fh = new BamUnAln(path.c_str());
    else if (in_file_type == FileT::bustard)
        fh = new Bustard(path.c_str());

    //
    // Open a file for recording discarded reads
    //
    if (discards) {
        path = path + ".discards";
        discard_fh = new ofstream(path.c_str(), ifstream::out);

        if (discard_fh->fail()) {
            cerr << "Error opening discard output file '" << path << "'\n";
            exit(1);
        }
    }

    //
    // Read in the first record, initializing the Seq object s. Then
    // initialize the Read object r, then loop, using the same objects.
    //
    Seq *s = fh->next_seq();
    if (s == NULL) {
        cerr << "Attempting to read first input record, unable to allocate "
             << "Seq object (Was the correct input type specified?).\n";
        exit(1);
    }

    r = new RawRead(strlen(s->seq), 1, min_bc_size_1, win_size);

    BarcodePair bc;
    //
    // If no barcodes were specified, set the barcode object to be the input file name so
    // that reads are written to an output file of the same name as the input file.
    //
    if (max_bc_size_1 == 0)
        bc.set(prefix);

    //cerr << "Length: " << r->len << "; Window length: " << r->win_len << "; Stop position: " << r->stop_pos << "\n";

    long i = 1;

    do {
        if (i % 10000 == 0) cerr << "  Processing short read " << i << "       \r";
        counter["total"]++;

        parse_input_record(s, r);

        if (barcode_type == inline_null ||
            barcode_type == index_null)
            bc.set(r->se_bc);
        else if (barcode_type == index_inline ||
                 barcode_type == inline_index)
            bc.set(r->se_bc, r->pe_bc);

        process_barcode(r, NULL, bc, pair_1_fhs, se_bc, pe_bc, barcode_log, counter);

        //
        // Adjust the size of the read to accommodate truncating the sequence and variable
        // barcode lengths. With standard Illumina data we want to output constant length
        // reads even as the barcode size may change. Other technologies, like IonTorrent
        // need to be truncated uniformly.
        //
        if (truncate_seq > 0) {
            if (truncate_seq + r->inline_bc_len <= r->len)
                r->set_len(truncate_seq + r->inline_bc_len);
        } else {
            if (barcode_type == inline_null || barcode_type == inline_inline ||        barcode_type == inline_index)
                r->set_len(r->len - (max_bc_size_1 - r->inline_bc_len));
        }

        if (r->retain)
            process_singlet(r, false, barcode_log[bc], counter);

        int result = 1;

        if (r->retain) {
            if (retain_header)
                result = (out_file_type == FileT::fastq || out_file_type == FileT::gzfastq) ?
                    write_fastq(pair_1_fhs[bc], s, r) :
                    write_fasta(pair_1_fhs[bc], s, r);
            else
                result = (out_file_type == FileT::fastq || out_file_type == FileT::gzfastq) ?
                    write_fastq(pair_1_fhs[bc], r, overhang) :
                    write_fasta(pair_1_fhs[bc], r, overhang);
        }

        if (!result) {
            cerr << "Error writing to output file for '" << bc.str() << "'\n";
            return_val = -1;
            break;
        }

        if (discards && !r->retain)
            result = (out_file_type == FileT::fastq || out_file_type == FileT::gzfastq) ?
                write_fastq(discard_fh, s) :
                write_fasta(discard_fh, s);

        if (!result) {
            cerr << "Error writing to discard file for '" << bc.str() << "'\n";
            return_val = -1;
            break;
        }

        delete s;
        i++;
    } while ((s = fh->next_seq()) != NULL);

    if (discards) delete discard_fh;

    //
    // Close the file and delete the Input object.
    //
    delete fh;

    return return_val;
}

inline int
process_singlet(RawRead *href,
                bool paired_end,
                map<string, long> &bc_log, map<string, long> &counter)
{
    if (filter_illumina && href->filter) {
        counter["ill_filtered"]++;
        href->retain = 0;
        return 0;
    }

    //
    // Drop this sequence if it has any uncalled nucleotides
    //
    if (clean) {
        for (char *p = href->seq + href->inline_bc_len; *p != '\0'; p++)
            if (*p == '.' || *p == 'N') {
                counter["low_quality"]++;
                href->retain = 0;
                return 0;
            }
    }

    bool adapter_trim = false;
    bool quality_trim = false;

    //
    // Drop or trim this sequence if it has low quality scores
    //
    if (quality) {
        int res = check_quality_scores(href, qual_offset, score_limit, len_limit, href->inline_bc_len);

        if (trim_reads) {
            if (res == 0) {
                counter["low_quality"]++;
                href->retain = 0;
                return 0;
            } else if (res < 0) {
                quality_trim = true;
            }
        } else {
            if (res <= 0) {
                counter["low_quality"]++;
                href->retain = 0;
                return 0;
            }
        }
    }

    //
    // Drop or trim this sequence if it contains adapter sequence.
    //
    if (filter_adapter) {
        int res = 1;
        if (paired_end == true  && adp_2_len > 0)
            res = filter_adapter_seq(href, adapter_2, adp_2_len, adp_2_kmers,
                                     kmer_size, distance, len_limit);
        if (paired_end == false && adp_1_len > 0)
            res = filter_adapter_seq(href, adapter_1, adp_1_len, adp_1_kmers,
                                     kmer_size, distance, len_limit);
        if (res == 0) {
            counter["adapter"]++;
            href->retain = 0;
            return 0;

        } else if (res < 0) {
            counter["adapter"]++;
            adapter_trim = true;
        }
    }

    if (adapter_trim || quality_trim)
        counter["trimmed"]++;

    if (barcode_type != null_null)
        bc_log["retained"]++;
    counter["retained"]++;

    return 0;
}

int dist(const char *res_enz, char *seq) {
    const char *p; char *q;

    int dist = 0;

    for (p = res_enz, q = seq; *p != '\0'; p++, q++)
        if (*p != *q) dist++;

    return dist;
}

int
print_results(int argc, char **argv,
              vector<BarcodePair> &barcodes,
              map<string, map<string, long> > &counters,
              map<BarcodePair, map<string, long> > &barcode_log)
{
    map<string, map<string, long> >::iterator it;

    string log_path = out_path + "process_shortreads.log";
    ofstream log(log_path.c_str());

    if (log.fail()) {
        cerr << "Unable to open log file '" << log_path << "'\n";
        return 0;
    }

    cerr << "Outputing details to log: '" << log_path << "'\n\n";

    init_log(log, argc, argv);

    log << "File\t"
        << "Retained Reads\t";
    if (filter_illumina)
        log << "Illumina Filtered\t";
    if (filter_adapter)
        log << "Adapter Seq" << "\t";
    log << "Low Quality\t"
        << "Ambiguous Barcodes\t"
        << "Trimmed Reads\t"
        << "Orphaned paired-end reads\t"
        << "Total\n";

    for (it = counters.begin(); it != counters.end(); it++) {
        log << it->first                 << "\t"
            << it->second["retained"]    << "\t";
        if (filter_illumina)
            log << it->second["ill_filtered"] << "\t";
        if (filter_adapter)
            log << it->second["adapter"] << "\t";
        log << it->second["low_quality"] << "\t"
            << it->second["ambiguous"]   << "\t"
            << it->second["trimmed"]    << "\t"
            << it->second["orphaned"]    << "\t"
            << it->second["total"]       << "\n";
    }

    map<string, long> c;
    c["total"]        = 0;
    c["low_quality"]  = 0;
    c["adapter"]      = 0;
    c["ill_filtered"] = 0;
    c["ambiguous"]    = 0;
    c["trimmed"]      = 0;
    c["orphaned"]     = 0;

    //
    // Total up the individual counters
    //
    for (it = counters.begin(); it != counters.end(); it++) {
        c["total"]        += it->second["total"];
        c["ill_filtered"] += it->second["ill_filtered"];
        c["adapter"]      += it->second["adapter"];
        c["low_quality"]  += it->second["low_quality"];
        c["ambiguous"]    += it->second["ambiguous"];
        c["trimmed"]      += it->second["trimmed"];
        c["orphaned"]     += it->second["orphaned"];
        c["retained"]     += it->second["retained"];
    }

    cerr << c["total"] << " total sequences;\n";
    if (filter_illumina)
        cerr << "  " << c["ill_filtered"] << " failed Illumina filtered reads;\n";
    if (filter_adapter)
        cerr << "  " << c["adapter"] << " reads contained adapter sequence;\n";
    cerr << "  " << c["ambiguous"]   << " ambiguous barcode drops;\n"
         << "  " << c["low_quality"] << " low quality read drops;\n"
         << "  " << c["trimmed"]     << " trimmed reads;\n"
         << "  " << c["orphaned"]    << " orphaned paired-end reads;\n"
         << c["retained"] << " retained reads.\n";

    log        << "\n"
        << "Total Sequences\t"      << c["total"]       << "\n";
    if (filter_illumina)
        log << "Failed Illumina filtered reads\t" << c["ill_filtered"] << "\n";
    if (filter_adapter)
        log << "Reads containing adapter sequence\t" << c["adapter"] << "\n";
    log
        << "Ambiguous Barcodes\t"   << c["ambiguous"]   << "\n"
        << "Low Quality\t"          << c["low_quality"] << "\n"
        << "Trimmed Reads\t"        << c["trimmed"]     << "\n"
        << "Orphaned Paired-ends\t" << c["orphaned"]    << "\n"
        << "Retained Reads\t"       << c["retained"]      << "\n";

    if (max_bc_size_1 == 0) return 0;

    //
    // Where barcode filenames specified?
    //
    bool bc_names = false;
    for (uint i = 0; i < barcodes.size(); i++)
        if (barcodes[i].name_exists()) {
            bc_names = true;
            break;
        }

    //
    // Print out barcode information.
    //
    log << "\n"
        << "Barcode\t";
    if (bc_names)
        log << "Filename\t";
    log << "Total\t"
        << "Retained\n";

    set<BarcodePair> barcode_list;

    for (uint i = 0; i < barcodes.size(); i++) {
        barcode_list.insert(barcodes[i]);

        log << barcodes[i] << "\t";
        if (bc_names)
            log << barcodes[i].name << "\t";
        if (barcode_log.count(barcodes[i]) == 0)
            log << "0\t" << "0\t" << "0\n";
        else
            log << barcode_log[barcodes[i]]["total"]    << "\t"
                << barcode_log[barcodes[i]]["retained"] << "\n";
    }

    log << "\n"
        << "Sequences not recorded\n"
        << "Barcode\t"
        << "Total\n";

    //
    // Sort unused barcodes by number of occurances.
    //
    map<BarcodePair, map<string, long> >::iterator bit;
    vector<pair<BarcodePair, int> > bcs;
    for (bit = barcode_log.begin(); bit != barcode_log.end(); bit++)
        bcs.push_back(make_pair(bit->first, bit->second["total"]));
    sort(bcs.begin(), bcs.end(), compare_barcodes);

    for (uint i = 0; i < bcs.size(); i++) {
        if (barcode_list.count(bcs[i].first)) continue;
        if (bcs[i].second == 0) continue;

        log << bcs[i].first << "\t"
            << bcs[i].second << "\n";
    }

    log.close();

    return 0;
}

int  compare_barcodes(pair<BarcodePair, int> a, pair<BarcodePair, int> b) {
    return a.second > b.second;
}

int parse_command_line(int argc, char* argv[]) {
    int c;

    while (1) {
        static struct option long_options[] = {
            {"help",                 no_argument, NULL, 'h'},
            {"version",              no_argument, NULL, 'v'},
            {"quality",              no_argument, NULL, 'q'},
            {"clean",                no_argument, NULL, 'c'},
            {"recover",              no_argument, NULL, 'r'},
            {"discards",             no_argument, NULL, 'D'},
            {"paired",               no_argument, NULL, 'P'},
            {"interleaved",          no_argument, NULL, 'I'},
            {"merge",                no_argument, NULL, 'm'},
            {"mate-pair",            no_argument, NULL, 'M'},
            {"no-overhang",          no_argument, NULL, 'O'}, {"no_overhang",          no_argument, NULL, 'O'},
            {"filter-illumina",      no_argument, NULL, 'F'}, {"filter_illumina",      no_argument, NULL, 'F'},
            {"retain-header",        no_argument, NULL, 'H'}, {"retain_header",        no_argument, NULL, 'H'},
            {"no-read-trimming",     no_argument, NULL, 'N'}, {"no_read_trimming",     no_argument, NULL, 'N'},
            {"null-index",           no_argument, NULL, 'U'}, {"null_index",           no_argument, NULL, 'U'},
            {"index-null",           no_argument, NULL, 'u'}, {"index_null",           no_argument, NULL, 'u'},
            {"inline-null",          no_argument, NULL, 'V'}, {"inline_null",          no_argument, NULL, 'V'},
            {"index-index",          no_argument, NULL, 'W'}, {"index_index",          no_argument, NULL, 'W'},
            {"inline-inline",        no_argument, NULL, 'x'}, {"inline_inline",        no_argument, NULL, 'x'},
            {"index-inline",         no_argument, NULL, 'Y'}, {"index_inline",         no_argument, NULL, 'Y'},
            {"inline-index",         no_argument, NULL, 'Z'}, {"inline_index",         no_argument, NULL, 'Z'},
            {"barcode-dist-1", required_argument, NULL, 'B'}, {"barcode_dist_1", required_argument, NULL, 'B'},
            {"barcode-dist-2", required_argument, NULL, 'C'}, {"barcode_dist_2", required_argument, NULL, 'C'},
            {"infile-type",    required_argument, NULL, 'i'}, {"infile_type",    required_argument, NULL, 'i'},
            {"outfile-type",   required_argument, NULL, 'y'}, {"outfile_type",   required_argument, NULL, 'y'},
            {"file",           required_argument, NULL, 'f'},
            {"file-p1",        required_argument, NULL, '1'}, {"file_p1",        required_argument, NULL, '1'},
            {"file-p2",        required_argument, NULL, '2'}, {"file_p2",        required_argument, NULL, '2'},
            {"path",           required_argument, NULL, 'p'},
            {"outpath",        required_argument, NULL, 'o'},
            {"truncate",       required_argument, NULL, 't'},
            {"barcodes",       required_argument, NULL, 'b'},
            {"window-size",    required_argument, NULL, 'w'}, {"window_size",    required_argument, NULL, 'w'},
            {"score-limit",    required_argument, NULL, 's'}, {"score_limit",    required_argument, NULL, 's'},
            {"encoding",       required_argument, NULL, 'E'},
            {"len-limit",      required_argument, NULL, 'L'}, {"len_limit",      required_argument, NULL, 'L'},
            {"adapter-1",      required_argument, NULL, 'A'}, {"adapter_1",      required_argument, NULL, 'A'},
            {"adapter-2",      required_argument, NULL, 'G'}, {"adapter_2",      required_argument, NULL, 'G'},
            {"adapter-mm",     required_argument, NULL, 'T'}, {"adapter_mm",     required_argument, NULL, 'T'},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;

        c = getopt_long(argc, argv, "hHvcqrINFuVWxYZOPmDi:y:f:o:t:B:C:b:1:2:p:s:w:E:L:A:G:T:", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            help();
            break;
        case 'i':
            if (strcasecmp(optarg, "bustard") == 0)
                in_file_type = FileT::bustard;
            else if (strcasecmp(optarg, "bam") == 0)
                in_file_type = FileT::bam;
            else if (strcasecmp(optarg, "gzfastq") == 0)
                in_file_type = FileT::gzfastq;
            else
                in_file_type = FileT::fastq;
            break;
        case 'y':
            if (strcasecmp(optarg, "fastq") == 0)
                out_file_type = FileT::fastq;
            else if (strcasecmp(optarg, "gzfastq") == 0)
                out_file_type = FileT::gzfastq;
            else if (strcasecmp(optarg, "fasta") == 0)
                out_file_type = FileT::fasta;
            else if (strcasecmp(optarg, "gzfasta") == 0)
                out_file_type = FileT::gzfasta;
            break;
        case 'E':
            if (strcasecmp(optarg, "phred64") == 0)
                qual_offset = 64;
            else if (strcasecmp(optarg, "phred33") == 0)
                qual_offset = 33;
            break;
        case 'f':
            in_file = optarg;
            break;
        case 'p':
            in_path_1 = optarg;
            in_path_2 = in_path_1;
            break;
        case '1':
            paired     = true;
            in_file_p1 = optarg;
            break;
        case '2':
            paired     = true;
            in_file_p2 = optarg;
            break;
        case 'P':
            paired = true;
            break;
        case 'I':
            interleaved = true;
            break;
        case 'B':
            barcode_dist_1 = is_integer(optarg);
            break;
        case 'C':
            barcode_dist_2 = is_integer(optarg);
            break;
        case 'o':
            out_path = optarg;
            break;
        case 'm':
            merge = true;
            break;
        case 'M':
            matepair = true;
            break;
        case 'D':
            discards = true;
            break;
        case 'q':
            quality = true;
            break;
        case 'c':
            clean = true;
            break;
        case 'r':
            recover = true;
            break;
        case 'O':
            overhang = false;
            break;
        case 'F':
            filter_illumina = true;
            break;
        case 'H':
            retain_header = true;
            break;
        case 'N':
            trim_reads = false;
            break;
        case 't':
            truncate_seq = is_integer(optarg);
            break;
        case 'b':
            barcode_file = optarg;
            if (barcode_type == null_null)
                barcode_type = inline_null;
            break;
        case 'U':
            barcode_type = null_index;
            break;
        case 'u':
            barcode_type = index_null;
            break;
        case 'V':
            barcode_type = inline_null;
            break;
        case 'W':
            barcode_type = index_index;
            break;
        case 'x':
            barcode_type = inline_inline;
            break;
        case 'Y':
            barcode_type = index_inline;
            break;
        case 'Z':
            barcode_type = inline_index;
            break;
        case 'A':
            adapter_1 = new char[strlen(optarg) + 1];
            strcpy(adapter_1, optarg);
            filter_adapter = true;
            break;
        case 'G':
            adapter_2 = new char[strlen(optarg) + 1];
            strcpy(adapter_2, optarg);
            filter_adapter = true;
            break;
        case 'T':
            distance = is_integer(optarg);
            break;
        case 'L':
            len_limit = is_integer(optarg);
            break;
        case 'w':
            win_size = is_double(optarg);
            break;
        case 's':
            score_limit = is_integer(optarg);
            break;
        case 'v':
            version();
            break;
        case '?':
            // getopt_long already printed an error message.
            help();
            break;

        default:
            cerr << "Unknown command line option '" << (char) c << "'\n";
            help();
            exit(1);
        }
    }

    if (optind < argc) {
        cerr << "Error: Failed to parse command line: '" << argv[optind] << "' is seen as a positional argument. Expected no positional arguments.\n";
        help();
    }

    if (in_file.length() == 0 && in_path_1.length() == 0 && in_file_p1.length() == 0) {
        cerr << "You must specify an input file of a directory path to a set of input files.\n";
        help();
    }

    if (in_file.length() > 0 && in_path_1.length() > 0) {
        cerr << "You must specify either a single input file (-f) or a directory path (-p), not both.\n";
        help();
    }

    if (in_file.length() > 0 && (in_file_p1.length() > 0 || in_file_p2.length() > 0)) {
        cerr << "You must specify either a single input file (-f) or a set of paired files (-1, -2), not both.\n";
        help();
    }

    if (in_path_1.length() > 0 && (in_file_p1.length() > 0 || in_file_p2.length() > 0)) {
        cerr << "You must specify either a file path (-p) or a set of paired files (-1, -2), not both.\n";
        help();
    }

    if (in_path_1.length() > 0 && in_path_1.at(in_path_1.length() - 1) != '/')
        in_path_1 += "/";

    if (in_path_2.length() > 0 && in_path_2.at(in_path_2.length() - 1) != '/')
        in_path_2 += "/";

    if (out_path.length() == 0)
        out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/')
        out_path += "/";

    if (barcode_file.length() == 0) {
        overhang = false;
        cerr << "No barcodes specified, files will not be demultiplexed.\n";
    }

    if (barcode_file.length() > 0 && merge) {
        cerr << "You may specify a set of barcodes, or that all files should be merged, not both.\n";
        help();
    }

    if (in_file_type == FileT::unknown)
        in_file_type = FileT::gzfastq;

    if (in_file_type == FileT::bam && paired == true && interleaved == false) {
        cerr << "You may only specify a BAM input file for paired-end data if the read pairs are interleaved.\n";
        help();
    }

    if (in_file_type == FileT::bam && (barcode_type != inline_null && barcode_type != inline_inline && barcode_type != null_null)) {
        cerr << "For BAM input files only inline or unbarcoded data can be processed.\n";
        help();
    }

    if (score_limit < 0 || score_limit > 40) {
        cerr << "Score limit must be between 0 and 40.\n";
        help();
    }

    if (win_size < 0 || win_size >= 1) {
        cerr << "Window size is a fraction between 0 and 1.\n";
        help();
    }

    if (recover && barcode_type != null_null) {
        if (barcode_type != index_null && barcode_type != inline_null && barcode_dist_2 < 0)
            barcode_dist_2 = barcode_dist_1;
    }

    return 0;
}

void version() {
    cerr << "process_shortreads " << VERSION << "\n\n";

    exit(1);
}

void help() {
    cerr << "process_shortreads " << VERSION << "\n"
              << "process_shortreads [-f in_file | -p in_dir [-P] [-I] | -1 pair_1 -2 pair_2] -b barcode_file -o out_dir [-i type] [-y type] [-c] [-q] [-r] [-E encoding] [-t len] [-D] [-w size] [-s lim] [-h]\n"
              << "  f: path to the input file if processing single-end seqeunces.\n"
              << "  i: input file type, either 'bustard' for the Illumina BUSTARD format, 'bam', 'fastq' (default), or 'gzfastq' for gzipped FASTQ.\n"
              << "  p: path to a directory of single-end Illumina files.\n"
              << "  1: first input file in a set of paired-end sequences.\n"
              << "  2: second input file in a set of paired-end sequences.\n"
              << "  P: specify that input is paired (for use with '-p').\n"
              << "  I: specify that the paired-end reads are interleaved in single files.\n"
              << "  o: path to output the processed files.\n"
              << "  y: output type, either 'fastq' or 'fasta' (default gzfastq).\n"
              << "  b: a list of barcodes for this run.\n"
              << "  c: clean data, remove any read with an uncalled base.\n"
              << "  q: discard reads with low quality scores.\n"
              << "  r: rescue barcodes.\n"
              << "  t: truncate final read length to this value.\n"
              << "  E: specify how quality scores are encoded, 'phred33' (Illumina 1.8+/Sanger, default) or 'phred64' (Illumina 1.3-1.5).\n"
              << "  D: capture discarded reads to a file.\n"
              << "  w: set the size of the sliding window as a fraction of the read length, between 0 and 1 (default 0.15).\n"
              << "  s: set the score limit. If the average score within the sliding window drops below this value, the read is discarded (default 10).\n"
              << "  h: display this help messsage.\n\n"
              << "  Barcode options:\n"
              << "    --inline-null:   barcode is inline with sequence, occurs only on single-end read (default).\n"
              << "    --index-null:    barcode is provded in FASTQ header (Illumina i5 or i7 read).\n"
              << "    --null-index:    barcode is provded in FASTQ header (Illumina i7 read if both i5 and i7 read are provided).\n"
              << "    --inline-inline: barcode is inline with sequence, occurs on single and paired-end read.\n"
              << "    --index-index:   barcode is provded in FASTQ header (Illumina i5 and i7 reads).\n"
              << "    --inline-index:  barcode is inline with sequence on single-end read and occurs in FASTQ header (from either i5 or i7 read).\n"
              << "    --index-inline:  barcode occurs in FASTQ header (Illumina i5 or i7 read) and is inline with single-end sequence (for single-end data) on paired-end read (for paired-end data).\n\n"
              << "  Adapter options:\n"
              << "    --adapter-1 <sequence>: provide adaptor sequence that may occur on the first read for filtering.\n"
              << "    --adapter-2 <sequence>: provide adaptor sequence that may occur on the paired-read for filtering.\n"
              << "      --adapter-mm <mismatches>: number of mismatches allowed in the adapter sequence.\n\n"
              << "  Output options:\n"
              << "    --retain-header: retain unmodified FASTQ headers in the output.\n"
              << "    --merge: if no barcodes are specified, merge all input files into a single output file (or single pair of files).\n\n"
              << "  Advanced options:\n"
              << "    --no-read-trimming: do not trim low quality reads, just discard them.\n"
              << "    --len-limit <limit>: when trimming sequences, specify the minimum length a sequence must be to keep it (default 31bp).\n"
              << "    --filter-illumina: discard reads that have been marked by Illumina's chastity/purity filter as failing.\n"
              << "    --barcode-dist: provide the distace between barcodes to allow for barcode rescue (default 2)\n"
              << "    --mate-pair: raw reads are circularized mate-pair data, first read will be reverse complemented.\n"
              << "    --no-overhang: data does not contain an overhang nucleotide between barcode and seqeunce.\n";

    exit(1);
}
