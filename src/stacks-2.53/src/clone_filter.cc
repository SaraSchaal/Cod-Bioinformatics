// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011-2017, Julian Catchen <jcatchen@illinois.edu>
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
// clone_filter -- find duplicate read pairs and reduce them to one representative
// pair of sequences in the data set. These reads are assumed to be the product of
// PCR amplification.
//

#include "clone_filter.h"

//
// Global variables to hold command-line options.
//
FileT    in_file_type  = FileT::unknown;
FileT    out_file_type = FileT::unknown;
string   in_file;
string   in_file_p1;
string   in_file_p2;
string   in_path_1;
string   in_path_2;
string   out_path;
bool     discards      = false;
bool     interleaved   = false;
bool     merge         = false;
bool     paired        = false;
bool     retain_oligo  = false;
barcodet barcode_type  = null_null;
int      oligo_len_1   = 0;
int      oligo_len_2   = 0;

//
// These variables are required for other linked objects, but we won't use them in clone_filter.
//
int    barcode_size;
uint   truncate_seq;
bool   ill_barcode;
bool   recover;
uint   min_bc_size_1;
uint   max_bc_size_1;
uint   min_bc_size_2;
uint   max_bc_size_2;
double win_size;

int main (int argc, char* argv[]) {
    IF_NDEBUG_TRY

    parse_command_line(argc, argv);

    if (oligo_len_2 == 0) oligo_len_2 = oligo_len_1;
    min_bc_size_1 = oligo_len_1;
    max_bc_size_1 = oligo_len_1;
    min_bc_size_2 = oligo_len_2;
    max_bc_size_2 = oligo_len_2;

    //
    // If input files are gzipped, output gziped files, unless the user chooses an output type.
    //
    if (out_file_type == FileT::unknown) {
        if (in_file_type == FileT::gzfastq || in_file_type == FileT::bam)
            out_file_type = FileT::gzfastq;
        else
            out_file_type = FileT::fastq;
    }

    if (paired)
        cerr << "Processing paired-end data.\n";
    else
        cerr << "Processing single-end data.\n";

    switch(barcode_type) {
    case null_null:
        cerr << "No oligo sequence specified, will use single and paired-end reads to determine clones.\n";
        break;
    case null_inline:
        cerr << "Searching for inline oligo on paired-end read.\n";
        break;
    case null_index:
        cerr << "Searching for index oligo (i7 Illumina read).\n";
        break;
    case index_null:
        cerr << "Searching for index oligo (i5 Illumina read).\n";
        break;
    case inline_null:
        cerr << "Searching for inline oligo on single-end read.\n";
        break;
    case index_index:
        cerr << "Searching for index oligos (i5 and i7 Illumina reads).\n";
        break;
    case inline_inline:
        cerr << "Searching for inline oligos on single and paired-end read.\n";
        break;
    case inline_index:
        cerr << "Searching for inline oligo on single-end read and index oligo (i5 or i7 Illumina read).\n";
        break;
    case index_inline:
        if (paired)
            cerr << "Searching for inline oligo on paired-end read and index oligo (i5 or i7 Illumina read).\n";
        else
            cerr << "Searching for inline oligo on single-end read and index oligo (i5 or i7 Illumina read).\n";
        break;
    }

    map<string, long> counters;
    counters["total"] = 0;
    counters["red_reads"] = 0;
    counters["dis_reads"] = 0;

    vector<pair<string, string> > files;

    build_file_list(files);

    CloneHash      clone_map;
    OligoHash      oligo_map;
    map<int, int>  clone_dist;
    vector<char *> clone_map_keys;

    int result = 1;
    for (uint i = 0; i < files.size(); i++) {
        cerr << "Processing file " << i+1 << " of " << files.size() << " [" << files[i].first.c_str() << "]\n";

        result = 1;
        if (paired) {
            if (barcode_type == null_null)
                result = process_paired_reads_by_sequence(files[i].first, files[i].second, counters, clone_map, clone_map_keys);
            else
                result = process_paired_reads(files[i].first, files[i].second, counters, oligo_map);

        } else {
            result = process_reads(files[i].first, counters, oligo_map);
        }

        if (!result) {
            cerr << "Error processing reads.\n";
            break;
        }
    }

    if (barcode_type == null_null && result) {
        write_clonereduced_sequence(files[0].first, files[0].second, clone_map, clone_dist, counters);

    } else {
        for (OligoHash::iterator i = oligo_map.begin(); i != oligo_map.end(); i++)
            for (map<string, uint16_t>::iterator j = i->second.begin(); j != i->second.end(); j++)
                clone_dist[j->second]++;
    }

    if (clone_map_keys.size() > 0) {
        cerr << "Freeing hash key memory...";
        free_hash(clone_map_keys);
        cerr << "done.\n";
    }

    //
    // Determine and print the distribution of read clones.
    //
    cerr << "Calculating the distribution of cloned read pairs...\n";

    vector<int> bins;
    map<int, int>::iterator it;
    for (it = clone_dist.begin(); it != clone_dist.end(); it++)
        bins.push_back(it->first);
    sort(bins.begin(), bins.end());
    cout << "Num Clones\tCount\n";
    for (uint i = 0; i < bins.size(); i++)
        cout << bins[i] << "\t" << clone_dist[bins[i]] << "\n";

    char buf[32];
    sprintf(buf, "%0.2f%%", ((double) (counters["total"] - counters["red_reads"]) / (double) counters["total"]) * 100);
    cerr << counters["total"] << " pairs of reads input. "
         << counters["red_reads"] << " pairs of reads output, discarded "
         << counters["dis_reads"] << " pairs of reads, " << buf << " clone reads.\n";

    return 0;
    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

int
process_paired_reads_by_sequence(string prefix_1, string prefix_2, map<string, long> &counters,
                                 CloneHash &clone_map, vector<char *> &clone_map_keys)
{
    Input *fh_1, *fh_2;

    int return_val = 1;

    string path_1 = in_path_1 + prefix_1;
    string path_2 = in_path_2 + prefix_2;

    cerr << "Reading data from:\n  "
         << path_1 << " and\n  "
         << path_2 << "\n";

    switch (in_file_type) {
    case FileT::fastq:
        fh_1 = new Fastq(path_1);
        fh_2 = interleaved ? fh_1 : new Fastq(path_2);
        break;
    case FileT::gzfastq:
        fh_1 = new GzFastq(path_1);
        fh_2 = interleaved ? fh_1 : new GzFastq(path_2);
        break;
    case FileT::fasta:
        fh_1 = new Fasta(path_1);
        fh_2 = interleaved ? fh_1 : new Fasta(path_2);
        break;
    case FileT::gzfasta:
        fh_1 = new GzFasta(path_1);
        fh_2 = interleaved ? fh_1 : new GzFasta(path_2);
        break;
    case FileT::bam:
        fh_1 = new BamUnAln(path_1);
        fh_2 = fh_1;
        break;
    case FileT::bustard:
        fh_1 = new Bustard(path_1);
        fh_2 = interleaved ? fh_1 : new Bustard(path_2);
    default:
        fh_1 = NULL;
        fh_2 = NULL;
        break;
    }

    //
    // Read in the first records, initializing the Seq objects, then loop, using the same objects.
    //
    Seq *s_1 = fh_1->next_seq();
    Seq *s_2 = fh_2->next_seq();
    if (s_1 == NULL || s_2 == NULL) {
        cerr << "Unable to allocate Seq object.\n";
        return 0;
    }

    long  i = 1;
    bool  exists;
    char *hash_key;
    uint  seq_len = strlen(s_1->seq);

    do {
        if (i % 10000 == 0) cerr << "Processing short read " << i << "       \r";

        counters["total"]++;

        exists = clone_map.count(s_1->seq) == 0 ? false : true;

        if (exists) {
            hash_key = s_1->seq;
        } else {
            hash_key = new char [seq_len + 1];
            strcpy(hash_key, s_1->seq);
            clone_map_keys.push_back(hash_key);
        }

        if (out_file_type == FileT::fastq ||
            out_file_type == FileT::gzfastq)
            clone_map[hash_key][s_2->seq].push_back(Pair(s_1->id, s_2->id, s_1->qual, s_2->qual));
        else if (out_file_type == FileT::fasta ||
                 out_file_type == FileT::gzfasta)
            clone_map[hash_key][s_2->seq].push_back(Pair(s_1->id, s_2->id));

        delete s_1;
        delete s_2;

        i++;
    } while ((s_1 = fh_1->next_seq()) != NULL &&
             (s_2 = fh_2->next_seq()) != NULL);

    cerr << "\n";

    delete fh_1;
    delete fh_2;

    return return_val;
}

int
write_clonereduced_sequence(string prefix_1, string prefix_2,
                            CloneHash &clone_map, map<int, int> &clone_dist,
                            map<string, long> &counters)
{
    ofstream  out_fh_1,   out_fh_2, discard_fh_1, discard_fh_2;
    gzFile    out_gzfh_1=NULL, out_gzfh_2=NULL, discard_gzfh_1=NULL, discard_gzfh_2=NULL;

    int return_val = 1;

    //
    // Open the input files.
    //
    string path_1;
    string path_2;

    //
    // Open the output files.
    //
    string suffix_1, suffix_2;

    if (out_file_type == FileT::gzfastq) {
        suffix_1 = ".1.fq.gz";
        suffix_2 = ".2.fq.gz";
    } else if (out_file_type == FileT::fastq) {
        suffix_1 = ".1.fq";
        suffix_2 = ".2.fq";
    } else if (out_file_type == FileT::gzfasta) {
        suffix_1 = ".1.fa.gz";
        suffix_2 = ".2.fa.gz";
    } else if (out_file_type == FileT::fasta) {
        suffix_1 = ".1.fa";
        suffix_2 = ".2.fa";
    }

    string file_1 = remove_suffix(in_file_type, prefix_1);
    path_1 = out_path + file_1 + suffix_1;
    if (out_file_type == FileT::gzfastq || out_file_type == FileT::gzfasta) {
        out_gzfh_1 = gzopen(path_1.c_str(), "wb");
        if (!(out_gzfh_1)) {
            cerr << "Error opening output file '" << path_1 << "'\n";
            return -1;
        }
    } else {
        out_fh_1.open(path_1.c_str(), ifstream::out);
        if (out_fh_1.fail()) {
            cerr << "Error opening output file '" << path_1 << "'\n";
            return -1;
        }
    }

    string file_2 = remove_suffix(in_file_type, prefix_2);
    path_2 = out_path + file_2 + suffix_2;
    if (out_file_type == FileT::gzfastq || out_file_type == FileT::gzfasta) {
        out_gzfh_2 = gzopen(path_2.c_str(), "wb");
        if (!(out_gzfh_2)) {
            cerr << "Error opening output file '" << path_2 << "'\n";
            return -1;
        }
    } else {
        out_fh_2.open(path_2.c_str(), ifstream::out);
        if (out_fh_2.fail()) {
            cerr << "Error opening output file '" << path_2 << "'\n";
            return -1;
        }
    }

    //
    // Open files for recording discarded reads.
    //
    if (discards) {
        path_1 = out_path + file_1 + ".discards" + suffix_1;
        if (out_file_type == FileT::gzfastq || out_file_type == FileT::gzfasta) {
            discard_gzfh_1 = gzopen(path_1.c_str(), "wb");
            if (!(discard_gzfh_1)) {
                cerr << "Error opening output file '" << path_1 << "'\n";
                return -1;
            }
        } else {
            discard_fh_1.open(path_1.c_str(), ifstream::out);
            if (discard_fh_1.fail()) {
                cerr << "Error opening discard output file '" << path_1 << "'\n";
                return -1;
            }
        }

        path_2 = out_path + file_2 + ".discards" + suffix_2;
        if (out_file_type == FileT::gzfastq || out_file_type == FileT::gzfasta) {
            discard_gzfh_2 = gzopen(path_2.c_str(), "wb");
            if (!(discard_gzfh_2)) {
                cerr << "Error opening output file '" << path_2 << "'\n";
                return -1;
            }
        } else {
            discard_fh_2.open(path_2.c_str(), ifstream::out);
            if (discard_fh_2.fail()) {
                cerr << "Error opening discard output file '" << path_2 << "'\n";
                return -1;
            }
        }
    }

    CloneHash::iterator hash_it;
    map<string, vector<Pair> >::iterator map_it;
    stringstream sstr_1, sstr_2;

    cerr << "Writing filtered data...\n";

    for (hash_it = clone_map.begin(); hash_it != clone_map.end(); hash_it++) {

        for (map_it = hash_it->second.begin(); map_it != hash_it->second.end(); map_it++) {
            sstr_1.str("");
            sstr_2.str("");

            if (out_file_type == FileT::gzfastq || out_file_type == FileT::fastq) {
                sstr_1 << "@" << map_it->second[0].p1_id << "\n"
                       << hash_it->first << "\n"
                       << "+\n"
                       << map_it->second[0].p1_qual << "\n";
                sstr_2 << "@" << map_it->second[0].p2_id << "\n"
                       << map_it->first << "\n"
                       << "+\n"
                       << map_it->second[0].p2_qual << "\n";
            } else {
                sstr_1 << ">" << map_it->second[0].p1_id << "\n"
                       << hash_it->first << "\n";
                sstr_2 << ">" << map_it->second[0].p2_id << "\n"
                       << map_it->first << "\n";
            }

            switch(out_file_type) {
            case FileT::gzfastq:
            case FileT::gzfasta:
                gzputs(out_gzfh_1, sstr_1.str().c_str());
                gzputs(out_gzfh_2, sstr_2.str().c_str());
                break;
            case FileT::fastq:
            case FileT::fasta:
            default:
                out_fh_1 << sstr_1.str();
                out_fh_2 << sstr_2.str();
            }

            counters["dis_reads"] += map_it->second.size() - 1;
            clone_dist[map_it->second.size()]++;

            //
            // Write cloned read pairs that we are discarding
            //
            if (discards) {
                for (uint i = 1; i < map_it->second.size(); i++) {
                    sstr_1.str("");
                    sstr_2.str("");

                    if (out_file_type == FileT::gzfastq || out_file_type == FileT::fastq) {
                        sstr_1 << "@" << map_it->second[i].p1_id << "\n"
                               << hash_it->first << "\n"
                               << "+\n"
                               << map_it->second[i].p1_qual << "\n";
                        sstr_2 << "@" << map_it->second[i].p2_id << "\n"
                               << map_it->first << "\n"
                               << "+\n"
                               << map_it->second[i].p2_qual << "\n";
                    } else {
                        sstr_1 << ">" << map_it->second[i].p1_id << "\n"
                               << hash_it->first << "\n";
                        sstr_2 << ">" << map_it->second[i].p2_id << "\n"
                               << map_it->first << "\n";
                    }

                    switch(out_file_type) {
                    case FileT::gzfastq:
                    case FileT::gzfasta:
                        gzputs(discard_gzfh_1, sstr_1.str().c_str());
                        gzputs(discard_gzfh_2, sstr_2.str().c_str());
                        break;
                    case FileT::fastq:
                    case FileT::fasta:
                    default:
                        discard_fh_1 << sstr_1.str();
                        discard_fh_2 << sstr_2.str();
                    }
                }
            }
            counters["red_reads"]++;
        }
    }

    cerr << "done.\n";

    if (out_file_type == FileT::gzfastq || out_file_type == FileT::gzfasta) {
        gzclose(out_gzfh_1);
        gzclose(out_gzfh_2);
        if (discards) {
            gzclose(discard_gzfh_1);
            gzclose(discard_gzfh_2);
        }
    } else {
        out_fh_1.close();
        out_fh_2.close();
        if (discards) {
            discard_fh_1.close();
            discard_fh_2.close();
        }
    }

    return return_val;
}
int
process_paired_reads(string prefix_1, string prefix_2, map<string, long> &counters, OligoHash &oligo_map)
{
    Input    *fh_1, *fh_2;
    RawRead  *r_1,  *r_2;
    ofstream  out_fh_1,   out_fh_2, discard_fh_1, discard_fh_2;
    gzFile    out_gzfh_1, out_gzfh_2, discard_gzfh_1, discard_gzfh_2;

    int return_val = 1;

    //
    // Open the input files.
    //
    string path_1 = in_path_1 + prefix_1;
    string path_2 = in_path_2 + prefix_2;

    if (interleaved)
        cerr << "  Reading data from:\n  " << path_1 << "\n";
    else
        cerr << "  Reading data from:\n  " << path_1 << " and\n  " << path_2 << "\n";

    switch (in_file_type) {
    case FileT::fastq:
        fh_1 = new Fastq(path_1);
        fh_2 = interleaved ? fh_1 : new Fastq(path_2);
        break;
    case FileT::gzfastq:
        fh_1 = new GzFastq(path_1);
        fh_2 = interleaved ? fh_1 : new GzFastq(path_2);
        break;
    case FileT::fasta:
        fh_1 = new Fasta(path_1);
        fh_2 = interleaved ? fh_1 : new Fasta(path_2);
        break;
    case FileT::gzfasta:
        fh_1 = new GzFasta(path_1);
        fh_2 = interleaved ? fh_1 : new GzFasta(path_2);
        break;
    case FileT::bam:
        fh_1 = new BamUnAln(path_1);
        fh_2 = fh_1;
        break;
    case FileT::bustard:
        fh_1 = new Bustard(path_1);
        fh_2 = interleaved ? fh_1 : new Bustard(path_2);
    default:
        fh_1 = NULL;
        fh_2 = NULL;
        break;
    }

    //
    // Open the output files.
    //
    string suffix_1, suffix_2;

    if (out_file_type == FileT::gzfastq) {
        suffix_1 = ".1.fq.gz";
        suffix_2 = ".2.fq.gz";
    } else if (out_file_type == FileT::fastq) {
        suffix_1 = ".1.fq";
        suffix_2 = ".2.fq";
    } else if (out_file_type == FileT::gzfasta) {
        suffix_1 = ".1.fa.gz";
        suffix_2 = ".2.fa.gz";
    } else if (out_file_type == FileT::fasta) {
        suffix_1 = ".1.fa";
        suffix_2 = ".2.fa";
    }

    string file_1 = remove_suffix(in_file_type, prefix_1);
    path_1 = out_path + file_1 + suffix_1;
    if (in_file_type == FileT::gzfastq || in_file_type == FileT::gzfasta) {
        out_gzfh_1 = gzopen(path_1.c_str(), "wb");
        if (!(out_gzfh_1)) {
            cerr << "Error opening output file '" << path_1 << "'\n";
            return -1;
        }
    } else {
        out_fh_1.open(path_1.c_str(), ifstream::out);
        if (out_fh_1.fail()) {
            cerr << "Error opening output file '" << path_1 << "'\n";
            return -1;
        }
    }

    string file_2 = remove_suffix(in_file_type, prefix_2);
    path_2 = out_path + file_2 + suffix_2;
    if (in_file_type == FileT::gzfastq || in_file_type == FileT::gzfasta) {
        out_gzfh_2 = gzopen(path_2.c_str(), "wb");
        if (!(out_gzfh_2)) {
            cerr << "Error opening output file '" << path_2 << "'\n";
            return -1;
        }
    } else {
        out_fh_2.open(path_2.c_str(), ifstream::out);
        if (out_fh_2.fail()) {
            cerr << "Error opening output file '" << path_2 << "'\n";
            return -1;
        }
    }

    //
    // Open files for recording discarded reads.
    //
    if (discards) {
        path_1 = out_path + file_1 + ".discards" + suffix_1;

        if (in_file_type == FileT::gzfastq || in_file_type == FileT::gzfasta) {
            discard_gzfh_1 = gzopen(path_1.c_str(), "wb");
            if (!(discard_gzfh_1)) {
                cerr << "Error opening discard file '" << path_1 << "'\n";
                return -1;
            }
        } else {
            discard_fh_1.open(path_1.c_str(), ifstream::out);
            if (discard_fh_1.fail()) {
                cerr << "Error opening discard file '" << path_1 << "'\n";
                return -1;
            }
        }

        path_2 = out_path + file_2 + ".discards" + suffix_2;

        if (in_file_type == FileT::gzfastq || in_file_type == FileT::gzfasta) {
            discard_gzfh_2 = gzopen(path_2.c_str(), "wb");
            if (!(discard_gzfh_2)) {
                cerr << "Error opening discard file '" << path_2 << "'\n";
                return -1;
            }
        } else {
            discard_fh_2.open(path_2.c_str(), ifstream::out);
            if (discard_fh_2.fail()) {
                cerr << "Error opening discard file '" << path_2 << "'\n";
                return -1;
            }
        }
    }

    //
    // Determine how much sequence we need to trim to remove the oligo seqeunce before printing.
    //
    int offset_1, offset_2;
    switch (barcode_type) {
    case null_inline:
        offset_1 = 0;
        offset_2 = oligo_len_1;
        break;
    case inline_null:
    case inline_index:
        offset_1 = oligo_len_1;
        offset_2 = 0;
        break;
    case null_index:
    case index_null:
    case index_index:
        offset_1 = 0;
        offset_2 = 0;
        break;
    case inline_inline:
        offset_1 = oligo_len_1;
        offset_2 = oligo_len_2;
        break;
    case index_inline:
        offset_1 = 0;
        offset_2 = oligo_len_2;
    default:
        DOES_NOT_HAPPEN;
        offset_1 = offset_2 = -1;
        break;
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

    long i        = 1;
    int  result_1 = 1;
    int  result_2 = 1;
    bool clone    = false;
    string oligo_1, oligo_2, key, oligo;

    do {
        if (i % 10000 == 0) cerr << "  Processing RAD-Tag " << i << "       \r";

        parse_input_record(s_1, r_1);
        parse_input_record(s_2, r_2);
        counters["total"]++;

        result_1 = 1;
        result_2 = 1;
        clone    = false;

        //
        // Fetch the randomized oligo sequence from the proper position in the reads.
        //
        switch (barcode_type) {
        case null_inline:
            oligo_1 = r_2->inline_bc;
            break;
        case inline_null:
            oligo_1 = r_1->inline_bc;
            break;
        case index_null:
            oligo_1 = r_1->index_bc;
            break;
        case null_index:
            oligo_1 = r_2->index_bc;
            break;
        case inline_inline:
            oligo_1 = r_1->inline_bc;
            oligo_2 = r_2->inline_bc;
            break;
        case index_index:
            oligo_1 = r_1->index_bc;
            oligo_2 = r_2->index_bc;
            break;
        case inline_index:
            oligo_1 = r_1->inline_bc;
            oligo_2 = r_2->index_bc;
            break;
        case index_inline:
            oligo_1 = r_1->index_bc;
            oligo_2 = r_2->inline_bc;
        default:
            break;
        }

        //
        // Have we seen this combination of oligos before for this read?
        //
        oligo = oligo_1 + oligo_2;
        key   = string(s_1->seq + offset_1) + string(s_2->seq + offset_2);

        // cerr << "Oligo: '" << oligo << "'\n"
        //      << "Seq: '" << s_1->seq << "'\n"
        //      << "Key: '" << key << "'\n";

        if (oligo_map.count(key) == 0)
            oligo_map[key] = map<string, uint16_t>();

        if (oligo_map[key].count(oligo) == 0) {
            oligo_map[key][oligo] = 1;
            clone = false;
        } else {
            oligo_map[key][oligo]++;
            clone = true;
        }

        if (clone == false) {
            counters["red_reads"]++;

            switch (out_file_type) {
            case FileT::fastq:
                result_1 = write_fastq(&out_fh_1, s_1, retain_oligo ? 0 : offset_1);
                result_2 = write_fastq(&out_fh_2, s_2, retain_oligo ? 0 : offset_2);
                break;
            case FileT::gzfastq:
                result_1 = write_fastq(&out_gzfh_1, s_1, retain_oligo ? 0 : offset_1);
                result_2 = write_fastq(&out_gzfh_2, s_2, retain_oligo ? 0 : offset_2);
                break;
            case FileT::fasta:
                result_1 = write_fasta(&out_fh_1, s_1, retain_oligo ? 0 : offset_1);
                result_2 = write_fasta(&out_fh_2, s_2, retain_oligo ? 0 : offset_2);
                break;
            case FileT::gzfasta:
                result_1 = write_fasta(&out_gzfh_1, s_1, retain_oligo ? 0 : offset_1);
                result_2 = write_fasta(&out_gzfh_2, s_2, retain_oligo ? 0 : offset_2);
            default:
                break;
            }

            if (!result_1 || !result_2) {
                cerr << "Error writing to output file for '" << file_1 << " / " << file_2 << "'\n";
                return_val = -1;
                break;
            }
        } else if (clone == true && discards) {
            counters["dis_reads"]++;

            switch (out_file_type) {
            case FileT::fastq:
                result_1 = write_fastq(&discard_fh_1,   s_1);
                result_2 = write_fastq(&discard_fh_2,   s_2);
                break;
            case FileT::gzfastq:
                result_1 = write_fastq(&discard_gzfh_1, s_1);
                result_2 = write_fastq(&discard_gzfh_2, s_2);
                break;
            case FileT::fasta:
                result_1 = write_fasta(&discard_fh_1,   s_1);
                result_2 = write_fasta(&discard_fh_2,   s_2);
                break;
            case FileT::gzfasta:
                result_1 = write_fasta(&discard_gzfh_1, s_1);
                result_2 = write_fasta(&discard_gzfh_2, s_2);
            default:
                break;
            }

            if (!result_1 || !result_2) {
                cerr << "Error writing to discard file for '" << file_1 << " / " << file_2 << "'\n";
                return_val = -1;
                break;
            }
        }

        delete s_1;
        delete s_2;

        i++;
    } while ((s_1 = fh_1->next_seq()) != NULL &&
             (s_2 = fh_2->next_seq()) != NULL);

    if (out_file_type == FileT::gzfastq || out_file_type == FileT::gzfasta) {
        gzclose(out_gzfh_1);
        gzclose(out_gzfh_2);
        if (discards) {
            gzclose(discard_gzfh_1);
            gzclose(discard_gzfh_2);
        }
    } else {
        out_fh_1.close();
        out_fh_2.close();
        if (discards) {
            discard_fh_1.close();
            discard_fh_2.close();
        }
    }

    delete fh_1;
    if (interleaved == false) delete fh_2;

    delete r_1;
    delete r_2;

    return return_val;
}

int
process_reads(string prefix_1, map<string, long> &counters, OligoHash &oligo_map)
{
    Input   *fh_1 = NULL;
    RawRead *r_1;
    ofstream out_fh_1, discard_fh_1;
    gzFile   out_gzfh_1, discard_gzfh_1;

    int return_val = 1;

    //
    // Open the input file.
    //
    string path_1 = in_path_1 + prefix_1;

    cerr << "  Reading data from:\n  " << path_1 << "\n";

    switch(in_file_type) {
    case FileT::fastq:
        fh_1 = new Fastq(path_1);
        break;
    case FileT::gzfastq:
        fh_1 = new GzFastq(path_1.c_str());
        break;
    case FileT::fasta:
        fh_1 = new Fasta(path_1);
        break;
    case FileT::gzfasta:
        fh_1 = new GzFasta(path_1);
        break;
    case FileT::bam:
        fh_1 = new BamUnAln(path_1);
        break;
    case FileT::bustard:
        fh_1 = new Bustard(path_1);
    default:
        break;
    }

    //
    // Open the output files.
    //
    string suffix_1;

    if (out_file_type == FileT::gzfastq)
        suffix_1 = ".fq.gz";
    else if (out_file_type == FileT::fastq)
        suffix_1 = ".fq";
    else if (out_file_type == FileT::gzfasta)
        suffix_1 = ".fa.gz";
    else if (out_file_type == FileT::fasta)
        suffix_1 = ".fa";

    string file_1 = prefix_1;
    int    pos    = file_1.find_last_of(".");
    if ((in_file_type == FileT::gzfastq || in_file_type == FileT::gzfasta) &&
        file_1.substr(pos) == ".gz") {
        file_1 = file_1.substr(0, pos);
        pos    = file_1.find_last_of(".");
    }
    path_1 = out_path + file_1.substr(0, pos) + suffix_1;
    if (in_file_type == FileT::gzfastq || in_file_type == FileT::gzfasta) {
        out_gzfh_1 = gzopen(path_1.c_str(), "wb");
        if (!(out_gzfh_1)) {
            cerr << "Error opening output file '" << path_1 << "'\n";
            return -1;
        }
    } else {
        out_fh_1.open(path_1.c_str(), ifstream::out);
        if (out_fh_1.fail()) {
            cerr << "Error opening output file '" << path_1 << "'\n";
            return -1;
        }
    }

    //
    // Open files for recording discarded reads.
    //
    if (discards) {
        path_1 = out_path + file_1 + ".discards" + suffix_1;

        if (in_file_type == FileT::gzfastq || in_file_type == FileT::gzfasta) {
            discard_gzfh_1 = gzopen(path_1.c_str(), "wb");
            if (!(discard_gzfh_1)) {
                cerr << "Error opening discard file '" << path_1 << "'\n";
                return -1;
            }
        } else {
            discard_fh_1.open(path_1.c_str(), ifstream::out);
            if (discard_fh_1.fail()) {
                cerr << "Error opening discard file '" << path_1 << "'\n";
                return -1;
            }
        }
    }

    //
    // Determine how much sequence we need to trim to remove the oligo seqeunce before printing.
    //
    int offset_1;
    switch (barcode_type) {
    case inline_null:
    case inline_index:
    case index_inline:
        offset_1 = oligo_len_1;
        break;
    default:
        offset_1 = 0;
        break;
    }

    //
    // Read in the first record, initializing the Seq object s. Then
    // initialize the Read object r, then loop, using the same objects.
    //
    Seq *s_1 = fh_1->next_seq();
    if (s_1 == NULL) {
        cerr << "Attempting to read first pair of input records, unable to allocate "
             << "Seq object (Was the correct input type specified?).\n";
        exit(1);
    }

    r_1 = new RawRead(strlen(s_1->seq), 1, min_bc_size_1, win_size);

    long   i        = 1;
    int    result_1 = 1;
    bool   clone    = false;
    string key, oligo_1;


    do {
        if (i % 10000 == 0) cerr << "  Processing RAD-Tag " << i << "       \r";

        parse_input_record(s_1, r_1);
        counters["total"]++;

        result_1 = 1;
        clone    = false;

        //
        // Fetch the randomized oligo sequence from the proper position in the reads.
        //
        if (barcode_type == inline_null)
            oligo_1 = r_1->inline_bc;
        else if (barcode_type == index_null)
            oligo_1 = r_1->index_bc;

        //
        // Have we seen this combination of oligos before for this read?
        //
        key = string(s_1->seq + offset_1);

        if (oligo_map.count(key) == 0)
            oligo_map[key] = map<string, uint16_t>();

        if (oligo_map[key].count(oligo_1) == 0) {
            oligo_map[key][oligo_1] = 1;
            clone = false;
        } else {
            oligo_map[key][oligo_1]++;
            clone = true;
        }

        if (clone == false) {
            counters["red_reads"]++;

            switch (out_file_type) {
            case FileT::fastq:
                result_1 = write_fastq(&out_fh_1,   s_1, retain_oligo ? 0 : offset_1);
                break;
            case FileT::gzfastq:
                result_1 = write_fastq(&out_gzfh_1, s_1, retain_oligo ? 0 : offset_1);
                break;
            case FileT::fasta:
                result_1 = write_fasta(&out_fh_1,   s_1, retain_oligo ? 0 : offset_1);
                break;
            case FileT::gzfasta:
                result_1 = write_fasta(&out_gzfh_1, s_1, retain_oligo ? 0 : offset_1);
            default:
                break;
            }

            if (!result_1) {
                cerr << "Error writing to output file for '" << file_1 << "'\n";
                return_val = -1;
                break;
            }
        } else if (clone == true && discards) {
            counters["dis_reads"]++;

            switch (out_file_type) {
            case FileT::fastq:
                result_1 = write_fastq(&discard_fh_1,   s_1);
                break;
            case FileT::gzfastq:
                result_1 = write_fastq(&discard_gzfh_1, s_1);
                break;
            case FileT::fasta:
                result_1 = write_fasta(&discard_fh_1,   s_1);
                break;
            case FileT::gzfasta:
                result_1 = write_fasta(&discard_gzfh_1, s_1);
            default:
                break;
            }

            if (!result_1) {
                cerr << "Error writing to discard file for '" << file_1 << "'\n";
                return_val = -1;
                break;
            }
        }

        delete s_1;

        i++;
    } while ((s_1 = fh_1->next_seq()) != NULL);

    if (out_file_type == FileT::gzfastq || out_file_type == FileT::gzfasta) {
        gzclose(out_gzfh_1);
        if (discards) gzclose(discard_gzfh_1);
    } else {
        out_fh_1.close();
        if (discards) discard_fh_1.close();
    }

    delete fh_1;

    delete r_1;

    return return_val;
}

int
free_hash(vector<char *> &keys)
{
    for (uint i = 0; i < keys.size(); i++) {
        delete [] keys[i];
    }
    keys.clear();

    return 0;
}

int parse_command_line(int argc, char* argv[]) {
    FileT ftype;
    int   c;

    while (1) {
        static struct option long_options[] = {
            {"help",          no_argument,       NULL, 'h'},
            {"version",       no_argument,       NULL, 'v'},
            {"discards",      no_argument,       NULL, 'D'},
            {"paired",        no_argument,       NULL, 'P'},
            {"null-index",    no_argument,       NULL, 'U'}, {"null_index",    no_argument,       NULL, 'U'},
            {"null-inline",   no_argument,       NULL, 'X'}, {"null_inline",   no_argument,       NULL, 'X'},
            {"index-null",    no_argument,       NULL, 'u'}, {"index_null",    no_argument,       NULL, 'u'},
            {"inline-null",   no_argument,       NULL, 'V'}, {"inline_null",   no_argument,       NULL, 'V'},
            {"index-index",   no_argument,       NULL, 'W'}, {"index_index",   no_argument,       NULL, 'W'},
            {"inline-inline", no_argument,       NULL, 'x'}, {"inline_inline", no_argument,       NULL, 'x'},
            {"index-inline",  no_argument,       NULL, 'Y'}, {"index_inline",  no_argument,       NULL, 'Y'},
            {"inline-index",  no_argument,       NULL, 'Z'}, {"inline_index",  no_argument,       NULL, 'Z'},
            {"infile-type",   required_argument, NULL, 'i'}, {"infile_type",   required_argument, NULL, 'i'},
            {"outfile-type",  required_argument, NULL, 'y'}, {"outfile_type",  required_argument, NULL, 'y'},
            {"file",          required_argument, NULL, 'f'},
            {"path",          required_argument, NULL, 'p'},
            {"file-p1",       required_argument, NULL, '1'}, {"file_p1",       required_argument, NULL, '1'},
            {"file-p2",       required_argument, NULL, '2'}, {"file_p2",       required_argument, NULL, '2'},
            {"outpath",       required_argument, NULL, 'o'},
            {"oligo-len-1",   required_argument, NULL, 'O'}, {"oligo_len_1",   required_argument, NULL, 'O'},
            {"oligo-len-2",   required_argument, NULL, 'L'}, {"oligo_len_2",   required_argument, NULL, 'L'},
            {"retain-oligo",  required_argument, NULL, 'R'}, {"retain_oligo",  required_argument, NULL, 'R'},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;

        c = getopt_long(argc, argv, "hvDPuUVWXxYZi:y:f:p:1:2:o:O:L:R:", long_options, &option_index);

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
            else if (strcasecmp(optarg, "fasta") == 0)
                in_file_type = FileT::fasta;
             else if (strcasecmp(optarg, "gzfasta") == 0)
                 in_file_type = FileT::gzfasta;
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
        case 'D':
            discards = true;
            break;
        case 'f':
            in_file = optarg;
            ftype   = FileT::fastq;
            break;
        case 'p':
            in_path_1 = optarg;
            in_path_2 = in_path_1;
            ftype     = FileT::fastq;
            break;
        case '1':
            paired     = true;
            in_file_p1 = optarg;
            ftype      = FileT::fastq;
            break;
        case '2':
            paired     = true;
            in_file_p2 = optarg;
            ftype      = FileT::fastq;
            break;
        case 'P':
            paired = true;
            break;
        case 'o':
            out_path = optarg;
            break;
        case 'U':
            barcode_type = null_index;
            break;
        case 'X':
            barcode_type = null_inline;
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
        case 'O':
            oligo_len_1 = is_integer(optarg);
            break;
        case 'L':
            oligo_len_2 = is_integer(optarg);
            break;
        case 'R':
            retain_oligo = true;
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

    if (in_file_type == FileT::unknown)
        in_file_type = ftype;

    if (paired == false && barcode_type == null_null) {
        cerr << "You must specify paired-end data if you do not have oligo sequences to differentiate cloned reads.\n";
        help();
    }

    if (barcode_type != null_null && oligo_len_1 == 0 && oligo_len_2 == 0) {
        cerr << "You must specify the length of the oligo sequences (--oligo-len-1 / --oligo-len-2).\n";
        help();
    }

    return 0;
}

void version() {
    std::cerr << "clone_filter " << VERSION << "\n\n";

    exit(1);
}

void help() {
    std::cerr << "clone_filter " << VERSION << "\n"
              << "clone_filter [-f in_file | -p in_dir [-P] [-I] | -1 pair_1 -2 pair_2] -o out_dir [-i type] [-y type] [-D] [-h]\n"
              << "  f: path to the input file if processing single-end sequences.\n"
              << "  p: path to a directory of files.\n"
              << "  P: files contained within directory specified by '-p' are paired.\n"
              << "  1: first input file in a set of paired-end sequences.\n"
              << "  2: second input file in a set of paired-end sequences.\n"
              << "  i: input file type, either 'bustard', 'fastq', 'fasta', 'gzfasta', or 'gzfastq' (default 'fastq').\n"
              << "  o: path to output the processed files.\n"
              << "  y: output type, either 'fastq', 'fasta', 'gzfasta', or 'gzfastq' (default same as input type).\n"
              << "  D: capture discarded reads to a file.\n"
              << "  h: display this help messsage.\n"
              << "  --oligo-len-1 len: length of the single-end oligo sequence in data set.\n"
              << "  --oligo-len-2 len: length of the paired-end oligo sequence in data set.\n"
              << "  --retain-oligo: do not trim off the random oligo sequence (if oligo is inline).\n\n"
              << "  Oligo sequence options:\n"
              << "    --inline-null:   random oligo is inline with sequence, occurs only on single-end read (default).\n"
              << "    --null-inline:   random oligo is inline with sequence, occurs only on the paired-end read.\n"
              << "    --null-index:    random oligo is provded in FASTQ header (Illumina i7 read if both i5 and i7 read are provided).\n"
              << "    --index-null:    random oligo is provded in FASTQ header (Illumina i5 or i7 read).\n"
              << "    --inline-inline: random oligo is inline with sequence, occurs on single and paired-end read.\n"
              << "    --index-index:   random oligo is provded in FASTQ header (Illumina i5 and i7 read).\n"
              << "    --inline-index:  random oligo is inline with sequence on single-end read and second oligo occurs in FASTQ header.\n"
              << "    --index-inline:  random oligo occurs in FASTQ header (Illumina i5 or i7 read) and is inline with sequence on single-end read (if single read data) or paired-end read (if paired data).\n\n";

    exit(1);
}
