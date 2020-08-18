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
// kmer_filter --
//

#include "kmer_filter.h"

//
// Global variables to hold command-line options.
//
FileT  in_file_type  = FileT::unknown;
FileT  out_file_type = FileT::fastq;
vector<string> in_files;
vector<string> in_pair_files;
string in_path;
string out_path;
string k_freq_path;
bool   filter_mare_k     = false;
bool   filter_abundant_k = false;
bool   normalize         = false;
bool   discards     = false;
bool   write_k_freq = false;
bool   read_k_freq  = false;
bool   kmer_distr   = false;
bool   ill_barcode  = false;
uint   truncate_seq = 0;
int    transition_lim = 3;
int    normalize_lim  = 0;
int    kmer_len       = 15;
int    max_k_freq     = 20000;
double max_k_pct      = 0.80;
double min_k_pct      = 0.15;
int    min_lim        = 0;
int    max_lim        = 0;
int    num_threads    = 1;

int main (int argc, char* argv[]) {
    IF_NDEBUG_TRY

    parse_command_line(argc, argv);

    if (min_lim == 0)
        min_lim = (int) round((double) kmer_len * 0.80);

    cerr << "Using a kmer size of " << kmer_len << "\n";
    if (filter_mare_k) {
        cerr << "Filtering out reads by identifying rare kmers: On.\n"
             << "  A kmer is considered rare when its coverage is at " << min_k_pct * 100 << "% or below the median kmer coverage for the read.\n"
             << "  A read is dropped when it contains " << min_lim << " or more rare kmers in a row.\n";
    } else
        cerr << "Filtering out reads by identifying rare kmers: Off.\n";

    if (filter_abundant_k) {
        cerr << "Filtering out reads by identifying abundant kmers: On.\n"
             << "  Kmer is considered abundant when it occurs " << max_k_freq << " or more times.\n";
        if (max_lim == 0)
            cerr << "  A read is dropped when it contains " << max_k_pct * 100 << "% or more abundant kmers.\n";
        else
            cerr << "  A read is dropped when it contains " << max_lim << " or more abundant kmers.\n";
    } else
        cerr << "Filtering out reads by identifying abundant kmers: Off.\n";

    if (normalize)
        cerr << "Normalizing read depth: On.\n"
             << "  Read depth limit: " << normalize_lim << "x\n";
    else
        cerr << "Normalizing read depth: Off.\n";

    vector<pair<string, string> >   files, pair_files;
    map<string, map<string, long> > counters;
    SeqKmerHash    kmers;
    vector<char *> kmers_keys;

    build_file_list(in_files, files);
    cerr << "Found " << files.size() << " input file(s).\n";

    build_file_list(in_pair_files, pair_files);
    cerr << "Found " << pair_files.size() << " paired input file(s).\n";

    if (filter_mare_k || filter_abundant_k || kmer_distr || write_k_freq) {
        cerr << "Generating kmer distribution...\n";

        if (read_k_freq)
            read_kmer_freq(k_freq_path, kmers, kmers_keys);
        else
            populate_kmers(pair_files, files, kmers, kmers_keys);

        // double kmer_med, kmer_mad;
        // calc_kmer_median(kmers, kmer_med, kmer_mad);
        // cerr << "Median kmer frequency: " << kmer_med << "; median absolute deviation: " << kmer_mad << "\n";

        if (kmer_distr) {
            generate_kmer_dist(kmers);
            if (write_k_freq == false)
                exit(1);
        }

        if (write_k_freq) {
            write_kmer_freq(k_freq_path, kmers);
            exit(1);
        }

        cerr << "Filtering reads by kmer frequency...\n";

        for (uint i = 0; i < pair_files.size(); i += 2) {
            cerr << "Processing paired file " << i+1 << " of " << (pair_files.size() / 2) << " [" << pair_files[i].second << "]\n";

            counters[pair_files[i].second]["total"]      = 0;
            counters[pair_files[i].second]["retained"]   = 0;
            counters[pair_files[i].second]["rare_k"]     = 0;
            counters[pair_files[i].second]["abundant_k"] = 0;

            process_paired_reads(pair_files[i].first,
                                 pair_files[i].second,
                                 pair_files[i+1].first,
                                 pair_files[i+1].second,
                                 kmers,
                                 counters[pair_files[i].second]);

            cerr << "  "
                 << counters[pair_files[i].second]["total"] << " total reads; "
                 << "-" << counters[pair_files[i].second]["rare_k"] << " rare k-mer reads; "
                 << "-" << counters[pair_files[i].second]["abundant_k"] << " abundant k-mer reads; "
                 << counters[pair_files[i].second]["retained"] << " retained reads.\n";
        }

        for (uint i = 0; i < files.size(); i++) {
            cerr << "Processing file " << i+1 << " of " << files.size() << " [" << files[i].second << "]\n";

            counters[files[i].second]["total"]      = 0;
            counters[files[i].second]["retained"]   = 0;
            counters[files[i].second]["rare_k"]     = 0;
            counters[files[i].second]["abundant_k"] = 0;

            process_reads(files[i].first,
                          files[i].second,
                          kmers,
                          counters[files[i].second]);

            cerr << "  "
                 << counters[files[i].second]["total"] << " total reads; "
                 << "-" << counters[files[i].second]["rare_k"] << " rare k-mer reads; "
                 << "-" << counters[files[i].second]["abundant_k"] << " abundant k-mer reads; "
                 << counters[files[i].second]["retained"] << " retained reads.\n";
        }

        free_kmer_hash(kmers, kmers_keys);

        print_results(counters);
    }

    if (normalize) {
        cerr << "Normalizing read depth...\n";

        //
        // Add the remainder files from the previous step to the queue.
        //
        if (filter_mare_k || filter_abundant_k) {
            string file;
            int    pos;
            for (uint i = 0; i < pair_files.size(); i += 2) {
                file  = pair_files[i].second;
                pos   = file.find_last_of(".");
                if (file.substr(pos - 2, 2) == ".1")
                    pos -= 2;
                file  = file.substr(0, pos) + ".rem.fil";
                file += out_file_type == FileT::fastq ? ".fq" : ".fa";
                cerr << "Adding remainder file generated in previous step to queue, '" << file << "\n";
                files.push_back(make_pair(pair_files[i].first, file));
            }
        }

        for (uint i = 0; i < pair_files.size(); i += 2) {
            cerr << "Processing paired files " << i+1 << " of " << (pair_files.size() / 2)
                 << " [" << pair_files[i].second << " / " << pair_files[i+1].second << "]\n";

            counters[pair_files[i].second]["total"]    = 0;
            counters[pair_files[i].second]["retained"] = 0;
            counters[pair_files[i].second]["overep"]   = 0;

            normalize_paired_reads(pair_files[i].first,
                                   pair_files[i].second,
                                   pair_files[i+1].first,
                                   pair_files[i+1].second,
                                   kmers, kmers_keys,
                                   counters[pair_files[i].second]);

            cerr << "  "
                 << counters[pair_files[i].second]["total"] << " total reads; "
                 << "-" << counters[pair_files[i].second]["overep"] << " over-represented reads; "
                 << counters[pair_files[i].second]["retained"] << " retained reads.\n";
        }

        for (uint i = 0; i < files.size(); i++) {
            cerr << "Processing file " << i+1 << " of " << files.size() << " [" << files[i].second << "]\n";

            counters[files[i].second]["total"]    = 0;
            counters[files[i].second]["retained"] = 0;
            counters[files[i].second]["overep"]   = 0;

            normalize_reads(files[i].first,
                            files[i].second,
                            kmers, kmers_keys,
                            counters[files[i].second]);

            cerr << "  "
                 << counters[files[i].second]["total"] << " total reads; "
                 << "-" << counters[files[i].second]["overep"] << " over-represented reads; "
                 << counters[files[i].second]["retained"] << " retained reads.\n";
        }

        free_kmer_hash(kmers, kmers_keys);
    }

    return 0;
    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

int process_paired_reads(string in_path_1,
                         string in_file_1,
                         string in_path_2,
                         string in_file_2,
                         SeqKmerHash &kmers,
                         map<string, long> &counter) {
    Input    *fh_1=NULL, *fh_2=NULL;
    ofstream *discard_fh_1=NULL, *discard_fh_2=NULL;
    int       pos;
    string    path;

    string path_1 = in_path_1 + in_file_1;
    string path_2 = in_path_2 + in_file_2;

    if (in_file_type == FileT::fastq) {
        fh_1 = new Fastq(path_1);
        fh_2 = new Fastq(path_2);
    } else if (in_file_type == FileT::fasta) {
        fh_1 = new Fasta(path_1);
        fh_2 = new Fasta(path_2);
    } else if (in_file_type == FileT::gzfasta) {
        fh_1 = new GzFasta(path_1 + ".gz");
        fh_2 = new GzFasta(path_2 + ".gz");
    } else if (in_file_type == FileT::gzfastq) {
        fh_1 = new GzFastq(path_1 + ".gz");
        fh_2 = new GzFastq(path_2 + ".gz");
    } else if (in_file_type == FileT::bustard) {
        fh_1 = new Bustard(path_1);
        fh_2 = new Bustard(path_2);
    }

    //
    // Open the output files.
    //
    pos  = in_file_1.find_last_of(".");
    path = out_path + in_file_1.substr(0, pos) + ".fil" + in_file_1.substr(pos);
    ofstream *ofh_1 = new ofstream(path.c_str(), ifstream::out);
    if (ofh_1->fail()) {
        cerr << "Error opening filtered output file '" << path << "'\n";
        exit(1);
    }

    pos  = in_file_2.find_last_of(".");
    path = out_path + in_file_2.substr(0, pos) + ".fil" + in_file_2.substr(pos);
    ofstream *ofh_2  = new ofstream(path.c_str(), ifstream::out);
    if (ofh_2->fail()) {
        cerr << "Error opening filtered paired output file '" << path << "'\n";
        exit(1);
    }

    pos  = in_file_2.find_last_of(".");
    //
    // Pull the ".2" suffix off the paired file name to make the remainder file name.
    //
    if (in_file_2.substr(pos - 2, 2) == ".2")
        pos -= 2;
    path  = out_path + in_file_2.substr(0, pos) + ".rem.fil";
    path += out_file_type == FileT::fastq ? ".fq" : ".fa";
    ofstream *rem_fh = new ofstream(path.c_str(), ifstream::out);
    if (rem_fh->fail()) {
        cerr << "Error opening filtered remainder output file '" << path << "'\n";
        exit(1);
    }

    //
    // Open a file for recording discarded reads
    //
    if (discards) {
        pos  = in_file_1.find_last_of(".");
        path = out_path + in_file_1.substr(0, pos) + ".discards" + in_file_1.substr(pos);
        discard_fh_1 = new ofstream(path.c_str(), ifstream::out);

        if (discard_fh_1->fail()) {
            cerr << "Error opening discard output file '" << path << "'\n";
            exit(1);
        }
        pos  = in_file_2.find_last_of(".");
        path = out_path + in_file_2.substr(0, pos) + ".discards" + in_file_2.substr(pos);
        discard_fh_2 = new ofstream(path.c_str(), ifstream::out);

        if (discard_fh_2->fail()) {
            cerr << "Error opening discard output file '" << path << "'\n";
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
        cerr << "Unable to allocate Seq object.\n";
        exit(1);
    }

    int   rare_k, abundant_k, num_kmers, max_kmer_lim;
    bool  retain_1, retain_2;
    char *kmer = new char[kmer_len + 1];

    long i = 1;
    do {
        if (i % 10000 == 0) cerr << "  Processing short read pair " << i << "       \r";
        counter["total"] += 2;
        stringstream msg_1, msg_2;

        retain_1       = true;
        retain_2       = true;
        num_kmers      = strlen(s_1->seq) - kmer_len + 1;
        max_kmer_lim   = max_lim == 0 ? (int) round((double) num_kmers * max_k_pct) : max_lim;

        //
        // Drop the first sequence if it has too many rare or abundant kmers.
        //
        kmer_lookup(kmers, s_1->seq, kmer, num_kmers, rare_k, abundant_k);

        if (filter_mare_k && rare_k > 0) {
            counter["rare_k"]++;
            retain_1 = false;
            msg_1 << "rare_k_" << rare_k;
        }

        if (retain_1 && filter_abundant_k && abundant_k > max_kmer_lim) {
            counter["abundant_k"]++;
            retain_1 = false;
            msg_1 << "abundant_k_" << abundant_k;
        }

        rare_k       = 0;
        abundant_k   = 0;
        num_kmers    = strlen(s_2->seq) - kmer_len + 1;
        max_kmer_lim = max_lim == 0 ? (int) round((double) num_kmers * max_k_pct) : max_lim;

        //
        // Drop the second sequence if it has too many rare or abundant kmers.
        //
        kmer_lookup(kmers, s_2->seq, kmer, num_kmers, rare_k, abundant_k);

        if (filter_mare_k && rare_k > 0) {
            counter["rare_k"]++;
            retain_2 = false;
            msg_2 << "rare_k_" << rare_k;
        }

        if (retain_2 && filter_abundant_k && abundant_k > max_kmer_lim) {
            counter["abundant_k"]++;
            retain_2 = false;
            msg_2 << "abundant_k_" << abundant_k;
        }

        if (retain_1 && retain_2) {
            counter["retained"] += 2;
            out_file_type == FileT::fastq ?
                 write_fastq(ofh_1, s_1) : write_fasta(ofh_1, s_1);
            out_file_type == FileT::fastq ?
                 write_fastq(ofh_2, s_2) : write_fasta(ofh_2, s_2);
        }

        if (retain_1 && !retain_2) {
            counter["retained"]++;
            out_file_type == FileT::fastq ?
                 write_fastq(rem_fh, s_1) : write_fasta(rem_fh, s_1);
        }

        if (!retain_1 && retain_2) {
            counter["retained"]++;
            out_file_type == FileT::fastq ?
                 write_fastq(rem_fh, s_2) : write_fasta(rem_fh, s_2);
        }

        if (discards && !retain_1)
            out_file_type == FileT::fastq ?
                write_fastq(discard_fh_1, s_1, msg_1.str()) : write_fasta(discard_fh_1, s_1, msg_1.str());
        if (discards && !retain_2)
            out_file_type == FileT::fastq ?
                write_fastq(discard_fh_2, s_2, msg_2.str()) : write_fasta(discard_fh_2, s_2, msg_2.str());

        delete s_1;
        delete s_2;

        i++;
    } while ((s_1 = fh_1->next_seq()) != NULL &&
             (s_2 = fh_2->next_seq()) != NULL);

    delete [] kmer;

    if (discards) {
        delete discard_fh_1;
        delete discard_fh_2;
    }

    //
    // Close the file and delete the Input object.
    //
    delete fh_1;
    delete fh_2;
    delete ofh_1;
    delete ofh_2;
    delete rem_fh;

    return 0;
}

int process_reads(string in_path,
                  string in_file,
                  SeqKmerHash &kmers,
                  map<string, long> &counter) {
    Input    *fh=NULL;
    ofstream *discard_fh=NULL;
    int       pos;

    string path = in_path + in_file;

    if (in_file_type == FileT::fastq)
        fh = new Fastq(path);
    else if (in_file_type == FileT::fasta)
        fh = new Fasta(path);
    else if (in_file_type == FileT::gzfastq)
        fh = new GzFastq(path + ".gz");
    else if (in_file_type == FileT::gzfasta)
        fh = new GzFasta(path + ".gz");
    else if (in_file_type == FileT::bustard)
        fh = new Bustard(path);

    //
    // Open the output file.
    //
    pos  = in_file.find_last_of(".");
    path = out_path + in_file.substr(0, pos) + ".fil" + in_file.substr(pos);
    ofstream *out_fh = new ofstream(path.c_str(), ifstream::out);

    if (out_fh->fail()) {
        cerr << "Error opening output file '" << path << "'\n";
        exit(1);
    }

    //
    // Open a file for recording discarded reads
    //
    if (discards) {
        pos  = in_file.find_last_of(".");
        path = out_path + in_file.substr(0, pos) + ".discards" + in_file.substr(pos);
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
        cerr << "Unable to allocate Seq object.\n";
        exit(1);
    }

    int   rare_k, abundant_k, num_kmers, max_kmer_lim;
    bool  retain;
    char *kmer  = new char[kmer_len + 1];
    long  i     = 1;

    do {
        if (i % 10000 == 0) cerr << "  Processing short read " << i << "       \r";
        counter["total"]++;
        stringstream msg;

        //
        // Drop this sequence if it has too many rare or abundant kmers.
        //
        retain       = true;
        num_kmers    = strlen(s->seq) - kmer_len + 1;
        max_kmer_lim = max_lim == 0 ? (int) round((double) num_kmers * max_k_pct) : max_lim;

        kmer_lookup(kmers, s->seq, kmer, num_kmers, rare_k, abundant_k);

        if (filter_mare_k && rare_k > 0) {
            counter["rare_k"]++;
            retain = false;
            msg << "rare_k_" << rare_k;
        }

        if (retain && filter_abundant_k && abundant_k > max_kmer_lim) {
            counter["abundant_k"]++;
            retain = false;
            msg << "abundant_k_" << abundant_k;
        }

        if (retain) {
            counter["retained"]++;
            out_file_type == FileT::fastq ?
                 write_fastq(out_fh, s) : write_fasta(out_fh, s);
        }

        if (discards && !retain)
            out_file_type == FileT::fastq ?
                write_fastq(discard_fh, s, msg.str()) : write_fasta(discard_fh, s, msg.str());

        delete s;

        i++;
    } while ((s = fh->next_seq()) != NULL);

    delete [] kmer;

    if (discards) delete discard_fh;

    //
    // Close the file and delete the Input object.
    //
    delete fh;
    delete out_fh;

    return 0;
}

int
normalize_paired_reads(string in_path_1,
                       string in_file_1,
                       string in_path_2,
                       string in_file_2,
                       SeqKmerHash &kmers, vector<char *> &kmer_keys,
                       map<string, long> &counter)
{
    Input    *fh_1=NULL, *fh_2=NULL;
    ofstream *discard_fh_1=NULL, *discard_fh_2=NULL;
    ofstream *ofh_1=NULL, *ofh_2=NULL, *rem_fh=NULL;
    string    path_1, path_2;
    int       pos;

    if (filter_abundant_k || filter_mare_k) {
        //
        // If we already filtered the data, open the files we created in the output
        // directory to normalize.
        //
        pos    = in_file_1.find_last_of(".");
        path_1 = out_path + in_file_1.substr(0, pos) + ".fil" + in_file_1.substr(pos);

        pos    = in_file_2.find_last_of(".");
        path_2 = out_path + in_file_2.substr(0, pos) + ".fil" + in_file_2.substr(pos);

        if (in_file_type == FileT::fastq) {
            fh_1 = new Fastq(path_1);
            fh_2 = new Fastq(path_2);
        } else if (in_file_type == FileT::gzfastq) {
            fh_1 = new Fastq(path_1);
            fh_2 = new Fastq(path_2);
        } else if (in_file_type == FileT::fasta) {
            fh_1 = new Fasta(path_1);
            fh_2 = new Fasta(path_2);
        } else if (in_file_type == FileT::gzfasta) {
            fh_1 = new Fasta(path_1);
            fh_2 = new Fasta(path_2);
        }
    } else {
        //
        // Otherwise, open unmodified files.
        //
        path_1 = in_path_1 + in_file_1;
        path_2 = in_path_2 + in_file_2;

        if (in_file_type == FileT::fastq) {
            fh_1 = new Fastq(path_1);
            fh_2 = new Fastq(path_2);
        } else if (in_file_type == FileT::gzfastq) {
            fh_1 = new GzFastq(path_1 + ".gz");
            fh_2 = new GzFastq(path_2 + ".gz");
        } else if (in_file_type == FileT::fasta) {
            fh_1 = new Fasta(path_1);
            fh_2 = new Fasta(path_2);
        } else if (in_file_type == FileT::gzfasta) {
            fh_1 = new GzFasta(path_1 + ".gz");
            fh_2 = new GzFasta(path_2 + ".gz");
        } else if (in_file_type == FileT::bustard) {
            fh_1 = new Bustard(path_1);
            fh_2 = new Bustard(path_2);
        }
    }

    //
    // Open the output files.
    //
    if (filter_abundant_k || filter_mare_k) {
        pos    = in_file_1.find_last_of(".");
        path_1 = out_path + in_file_1.substr(0, pos) + ".fil.norm" + in_file_1.substr(pos);
        ofh_1  = new ofstream(path_1.c_str(), ifstream::out);

        if (ofh_1->fail()) {
            cerr << "Error opening normalized output file '" << path_1 << "'\n";
            exit(1);
        }

        pos    = in_file_2.find_last_of(".");
        path_2 = out_path + in_file_2.substr(0, pos) + ".fil.norm" + in_file_2.substr(pos);
        ofh_2  = new ofstream(path_2.c_str(), ifstream::out);

        if (ofh_2->fail()) {
            cerr << "Error opening normalized paired output file '" << path_2 << "'\n";
            exit(1);
        }

        if (in_file_2.substr(pos - 2, 2) == ".2")
            pos -= 2;
        path_2  = out_path + in_file_2.substr(0, pos) + ".fil.norm.rem";
        path_2 += out_file_type == FileT::fastq ? ".fq" : ".fa";
        rem_fh  = new ofstream(path_2.c_str(), ifstream::out);

        if (rem_fh->fail()) {
            cerr << "Error opening normalized remainder output file '" << path_2 << "'\n";
            exit(1);
        }

    } else {
        pos    = in_file_1.find_last_of(".");
        path_1 = out_path + in_file_1.substr(0, pos) + ".norm" + in_file_1.substr(pos);
        ofh_1  = new ofstream(path_1.c_str(), ifstream::out);

        if (ofh_1->fail()) {
            cerr << "Error opening normalized output file '" << path_1 << "'\n";
            exit(1);
        }

        pos    = in_file_2.find_last_of(".");
        path_2 = out_path + in_file_2.substr(0, pos) + ".norm" + in_file_2.substr(pos);
        ofh_2  = new ofstream(path_2.c_str(), ifstream::out);

        if (ofh_2->fail()) {
            cerr << "Error opening normalized paired output file '" << path_2 << "'\n";
            exit(1);
        }

        if (in_file_2.substr(pos - 2, 2) == ".2")
            pos -= 2;
        path_2  = out_path + in_file_2.substr(0, pos) + ".norm.rem";
        path_2 += out_file_type == FileT::fastq ? ".fq" : ".fa";
        rem_fh  = new ofstream(path_2.c_str(), ifstream::out);

        if (rem_fh->fail()) {
            cerr << "Error opening normalized remainder output file '" << path_2 << "'\n";
            exit(1);
        }
    }

    //
    // Open a file for recording discarded reads
    //
    if (discards) {
        pos = in_file_1.find_last_of(".");
        if (filter_abundant_k || filter_mare_k)
            path_1 = out_path + in_file_1.substr(0, pos) + ".fil.discards" + in_file_1.substr(pos);
        else
            path_1 = out_path + in_file_1.substr(0, pos) + ".discards" + in_file_1.substr(pos);
        discard_fh_1 = new ofstream(path_1.c_str(), ifstream::out);

        if (discard_fh_1->fail()) {
            cerr << "Error opening discard output file '" << path_1 << "'\n";
            exit(1);
        }

        pos = in_file_2.find_last_of(".");
        if (filter_abundant_k || filter_mare_k)
            path_2 = out_path + in_file_2.substr(0, pos) + ".fil.discards" + in_file_2.substr(pos);
        else
            path_2 = out_path + in_file_2.substr(0, pos) + ".discards" + in_file_2.substr(pos);
        discard_fh_2 = new ofstream(path_2.c_str(), ifstream::out);

        if (discard_fh_2->fail()) {
            cerr << "Error opening discard output file '" << path_1 << "'\n";
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
        cerr << "Unable to allocate Seq object.\n";
        exit(1);
    }

    int   num_kmers;
    bool  retain_1, retain_2;
    char *kmer = new char[kmer_len + 1];

    long i = 1;
    do {
        if (i % 10000 == 0) cerr << "  Processing short read pair " << i << "       \r";
        counter["total"] += 2;

        retain_1  = true;
        retain_2  = true;
        num_kmers = strlen(s_1->seq) - kmer_len + 1;

        //
        // Drop the first sequence if it has too many rare or abundant kmers.
        //
        retain_1 = normalize_kmer_lookup(kmers, s_1->seq, kmer, num_kmers, kmer_keys);

        num_kmers = strlen(s_2->seq) - kmer_len + 1;

        //
        // Drop the second sequence if it has too many rare or abundant kmers.
        //
        retain_2 = normalize_kmer_lookup(kmers, s_2->seq, kmer, num_kmers, kmer_keys);

        if (retain_1 && retain_2) {
            counter["retained"] += 2;
            out_file_type == FileT::fastq ?
                 write_fastq(ofh_1, s_1) : write_fasta(ofh_1, s_1);
            out_file_type == FileT::fastq ?
                 write_fastq(ofh_2, s_2) : write_fasta(ofh_2, s_2);
        } else {
            counter["overep"] +=2;
        }

        if (retain_1 && !retain_2) {
            counter["retained"]++;
            counter["overep"]++;
            out_file_type == FileT::fastq ?
                 write_fastq(rem_fh, s_1) : write_fasta(rem_fh, s_1);
        }

        if (!retain_1 && retain_2) {
            counter["retained"]++;
            counter["overep"]++;
            out_file_type == FileT::fastq ?
                 write_fastq(rem_fh, s_2) : write_fasta(rem_fh, s_2);
        }

        if (discards && !retain_1)
            out_file_type == FileT::fastq ?
                write_fastq(discard_fh_1, s_1) : write_fasta(discard_fh_1, s_1);
        if (discards && !retain_2)
            out_file_type == FileT::fastq ?
                write_fastq(discard_fh_2, s_2) : write_fasta(discard_fh_2, s_2);

        delete s_1;
        delete s_2;

        i++;
    } while ((s_1 = fh_1->next_seq()) != NULL &&
             (s_2 = fh_2->next_seq()) != NULL);

    delete [] kmer;

    if (discards) {
        delete discard_fh_1;
        delete discard_fh_2;
    }

    //
    // Close the file and delete the Input object.
    //
    delete fh_1;
    delete fh_2;
    delete ofh_1;
    delete ofh_2;
    delete rem_fh;

    return 0;
}

int
normalize_reads(string in_path,
                string in_file,
                SeqKmerHash &kmers, vector<char *> &kmer_keys,
                map<string, long> &counter)
{
    Input    *fh=NULL;
    ofstream *discard_fh=NULL;
    string    path;

    int pos = in_file.find_last_of(".");

    if (filter_abundant_k || filter_mare_k) {
        if (in_file.substr(pos - 4, 4) == ".fil")
            path = out_path + in_file;
        else
            path = out_path + in_file.substr(0, pos) + ".fil" + in_file.substr(pos);

        if (in_file_type == FileT::fastq)
            fh = new Fastq(path);
        else if (in_file_type == FileT::gzfastq)
            fh = new Fastq(path);
        else if (in_file_type == FileT::fasta)
            fh = new Fasta(path);
        else if (in_file_type == FileT::gzfasta)
            fh = new Fasta(path);
        else if (in_file_type == FileT::bustard)
            fh = new Bustard(path);

    } else {
        path = in_path + in_file;

        if (in_file_type == FileT::fastq)
            fh = new Fastq(path);
        else if (in_file_type == FileT::gzfastq)
            fh = new GzFastq(path + ".gz");
        else if (in_file_type == FileT::fasta)
            fh = new Fasta(path);
        else if (in_file_type == FileT::gzfasta)
            fh = new GzFasta(path + ".gz");
        else if (in_file_type == FileT::bustard)
            fh = new Bustard(path);
    }

    //
    // Open the output file.
    //
    // if (filter_abundant_k || filter_mare_k) {
    //         path = out_path + in_file.substr(0, pos) + ".norm" + in_file.substr(pos);
    // } else {
    //         path = out_path + in_file.substr(0, pos) + ".norm" + in_file.substr(pos);
    // }
    path = out_path + in_file.substr(0, pos) + ".norm" + in_file.substr(pos);
    ofstream *out_fh = new ofstream(path.c_str(), ifstream::out);

    if (out_fh->fail()) {
        cerr << "Error opening normalized output file '" << path << "'\n";
        exit(1);
    }

    //
    // Open a file for recording discarded reads
    //
    if (discards) {
        if (filter_abundant_k || filter_mare_k)
            path = out_path + in_file.substr(0, pos) + ".fil.discards" + in_file.substr(pos);
        else
            path = out_path + in_file.substr(0, pos) + ".discards" + in_file.substr(pos);
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
        cerr << "Unable to allocate Seq object.\n";
        exit(1);
    }

    int   num_kmers;
    bool  retain;
    char *kmer  = new char[kmer_len + 1];
    long  i     = 1;

    do {
        if (i % 10000 == 0) cerr << "  Processing short read " << i << "       \r";
        counter["total"]++;

        //
        // Drop this sequence if it has too many rare or abundant kmers.
        //
        retain    = true;
        num_kmers = strlen(s->seq) - kmer_len + 1;

        retain = normalize_kmer_lookup(kmers, s->seq, kmer, num_kmers, kmer_keys);

        if (retain) {
            counter["retained"]++;
            out_file_type == FileT::fastq ?
                 write_fastq(out_fh, s) : write_fasta(out_fh, s);
        } else {
            counter["overep"]++;
        }

        if (discards && !retain)
            out_file_type == FileT::fastq ?
                write_fastq(discard_fh, s) : write_fasta(discard_fh, s);

        delete s;

        i++;
    } while ((s = fh->next_seq()) != NULL);

    delete [] kmer;

    if (discards) delete discard_fh;

    //
    // Close the file and delete the Input object.
    //
    delete fh;
    delete out_fh;

    return 0;
}

int
populate_kmers(vector<pair<string, string> > &pair_files,
               vector<pair<string, string> > &files,
               SeqKmerHash &kmers,
               vector<char *> &kmers_keys)
{
    //
    // Break each read down into k-mers and create a hash map of those k-mers
    // recording in which sequences they occur.
    //
    uint j   = 1;
    uint cnt = files.size() + pair_files.size();
    for (uint i = 0; i < files.size(); i++) {
        cerr << "Generating kmers from file " << j << " of " << cnt << " [" << files[i].second << "]\n";
        process_file_kmers(files[i].first + files[i].second, kmers, kmers_keys);
        j++;
    }

    for (uint i = 0; i < pair_files.size(); i++) {
        cerr << "Generating kmers from file " << j << " of " << cnt << " [" << pair_files[i].second << "]\n";
        process_file_kmers(pair_files[i].first + pair_files[i].second, kmers, kmers_keys);
        j++;
    }

    cerr << kmers.size() << " unique k-mers recorded.\n";

    return 0;
}

int
read_kmer_freq(string in_path, SeqKmerHash &kmer_map, vector<char *> &kmer_map_keys)
{
    cerr << "Reading kmer frequencies from '" << in_path.c_str() << "'...\n";

    ifstream fh(in_path.c_str(), ifstream::in);
    if (fh.fail()) {
        cerr << "Error opening rare kmer frequency input file '" << in_path << "'\n";
        exit(1);
    }

    char *hash_key;
    bool  exists;
    int   len, cnt;
    char  kmer[id_len];
    char  line[max_len];
    vector<string> parts;

    long i = 1;

    while (fh.good()) {
        if (i % 10000 == 0) cerr << "  Processing kmer " << i << "       \r";

        fh.getline(line, max_len);

        len = strlen(line);
        if (len == 0) continue;

        //
        // Check that there is no carraige return in the buffer.
        //
        if (line[len - 1] == '\r') line[len - 1] = '\0';

        //
        // Ignore comments
        //
        if (line[0] == '#') continue;

        //
        // Parse the kmer and the number of times it occurred
        // <kmer> <tab> <integer>
        //
        parse_tsv(line, parts);

        if (parts.size() != 2) {
            cerr << "kmer frequencies are not formated correctly: expecting two, tab separated columns, found " << parts.size() << ".\n";
            exit(1);
        }

        strcpy(kmer, parts[1].c_str());
        cnt = is_integer(kmer);
        if (cnt < 0) {
            cerr << "Non integer found in second column.\n";
            exit(1);
        }

        strcpy(kmer, parts[0].c_str());
        exists = kmer_map.count(kmer) == 0 ? false : true;

        if (exists) {
            cerr << "Warning: kmer '" << kmer << "' already exists in the kmer hash map.\n";
            hash_key = kmer;
            kmer_map[hash_key] += cnt;
        } else {
            hash_key = new char [strlen(kmer) + 1];
            strcpy(hash_key, kmer);
            kmer_map_keys.push_back(hash_key);
            kmer_map[hash_key]  = cnt;
        }
        i++;
    }

    fh.close();

    cerr << kmer_map.size() << " unique k-mers read.\n";

    kmer_len = strlen(kmer_map.begin()->first);
    cerr << "Setting kmer length to " << kmer_len << "bp.\n";

    return 0;
}

int
write_kmer_freq(string path, SeqKmerHash &kmer_map)
{
    cerr << "Writing kmer frequencies to '" << path.c_str() << "'...";

    ofstream out_fh(path.c_str(), ifstream::out);
    if (out_fh.fail()) {
        cerr << "Error opening rare kmer output file '" << path << "'\n";
        exit(1);
    }

    SeqKmerHash::iterator i;

    out_fh << "# Kmer\tCount\n";

    for (i = kmer_map.begin(); i != kmer_map.end(); i++) {
        out_fh << i->first << "\t" << i->second << "\n";
    }

    out_fh.close();

    cerr << "done.\n";

    return 0;
}

int
process_file_kmers(string path, SeqKmerHash &kmer_map, vector<char *> &kmer_map_keys)
{
    vector<char *> kmers;
    char          *hash_key;
    bool           exists;
    int            j;
    Input         *fh = NULL;

    if (in_file_type == FileT::fastq)
        fh = new Fastq(path);
    else if (in_file_type == FileT::gzfastq)
        fh = new GzFastq(path + ".gz");
    else if (in_file_type == FileT::fasta)
        fh = new Fasta(path);
    else if (in_file_type == FileT::gzfasta)
        fh = new GzFasta(path + ".gz");
    else if (in_file_type == FileT::bustard)
        fh = new Bustard(path.c_str());

    //
    // Read in the first record, initializing the Seq object s.
    //
    Seq *s = fh->next_seq();
    if (s == NULL) {
        cerr << "Unable to allocate Seq object.\n";
        exit(1);
    }

    int   num_kmers;
    char *kmer = new char [kmer_len + 1];

    long i = 1;
    do {
        if (i % 10000 == 0) cerr << "  Processing short read " << i << "       \r";

        num_kmers = strlen(s->seq) - kmer_len + 1;

        //
        // Generate and hash the kmers for this raw read
        //
        kmer[kmer_len] = '\0';

        for (j = 0; j < num_kmers; j++) {
            strncpy(kmer, s->seq + j, kmer_len);

            exists = kmer_map.count(kmer) == 0 ? false : true;

            if (exists) {
                    hash_key = kmer;
            } else {
                    hash_key = new char [kmer_len + 1];
                    strcpy(hash_key, kmer);
                kmer_map_keys.push_back(hash_key);
            }

            kmer_map[hash_key]++;
        }

        delete s;

        i++;
    } while ((s = fh->next_seq()) != NULL);

    delete [] kmer;

    //
    // Close the file and delete the Input object.
    //
    delete fh;

    return 0;
}

int
generate_kmer_dist(SeqKmerHash &kmer_map)
{
    SeqKmerHash::iterator i;
    map<uint, uint> bins;

    cerr << "Generating kmer distribution...\n";

    for (i = kmer_map.begin(); i != kmer_map.end(); i++)
        bins[i->second]++;

    map<uint, uint>::iterator j;
    vector<pair<uint, uint> > sorted_kmers;

    for (j = bins.begin(); j != bins.end(); j++)
        sorted_kmers.push_back(make_pair(j->first, j->second));

    cout << "KmerFrequency\tCount\n";

    for (unsigned long k = 0; k < sorted_kmers.size(); k++)
        cout << sorted_kmers[k].first << "\t" << sorted_kmers[k].second << "\n";

    return 0;
}

int
calc_kmer_median(SeqKmerHash &kmers, double &kmer_med, double &kmer_mad)
{
    kmer_med = 0.0;
    kmer_mad = 0.0;

    int num_kmers = kmers.size();
    vector<int> freqs, residuals;
    freqs.reserve(num_kmers);

    SeqKmerHash::iterator i;

    for (i = kmers.begin(); i != kmers.end(); i++)
        freqs.push_back(i->second);

    sort(freqs.begin(), freqs.end());

    kmer_med = num_kmers % 2 == 0 ?
        (double) (freqs[num_kmers / 2 - 1] + freqs[num_kmers / 2]) / 2.0 :
        (double) freqs[num_kmers / 2 - 1];

    //
    // Calculate the median absolute deviation.
    //
    residuals.reserve(num_kmers);

    for (int j = 0; j < num_kmers; j++)
        residuals.push_back(abs(freqs[j] - (int) kmer_med));
    sort(residuals.begin(), residuals.end());

    kmer_mad = num_kmers % 2 == 0 ?
        (double) (residuals[num_kmers / 2 - 1] + residuals[num_kmers / 2]) / 2.0 :
        (double) residuals[num_kmers / 2 - 1];

    return 0;
}

int
kmer_map_cmp(pair<char *, long> a, pair<char *, long> b)
{
    return (a.second < b.second);
}

inline bool
normalize_kmer_lookup(SeqKmerHash &kmer_map,
                      char *read, char *kmer,
                      int num_kmers,
                      vector<char *> &kmer_keys)
{
    kmer[kmer_len] = '\0';
    int  cnt       = 0;
    bool retain    = true;
    //
    // Generate kmers from this read, increment kmer frequency in dataset.
    //
    vector<int> sorted_cnts;
    sorted_cnts.reserve(num_kmers);

    // cout << "# " << read << "\n";

    for (int j = 0; j < num_kmers; j++) {
        strncpy(kmer, read + j, kmer_len);

        cnt = kmer_map.count(kmer) > 0 ? kmer_map[kmer] : 0;
        sorted_cnts.push_back(cnt);

        // cout << kmer << "\t" << j << "\t" << cnt << "\n";
    }

    //
    // Calculate the median kmer frequency along the read.
    //
    sort(sorted_cnts.begin(), sorted_cnts.end());
    double median = num_kmers % 2 == 0 ?
        (double) (sorted_cnts[num_kmers / 2 - 1] + sorted_cnts[num_kmers / 2]) / 2.0 :
        (double) sorted_cnts[num_kmers / 2 - 1];
    // cout << "# median: " << median << "\n";

    if (median > normalize_lim)
        retain = false;

    //
    // Generate and hash the kmers for this raw read
    //
    bool  exists;
    char *hash_key;
    kmer[kmer_len] = '\0';

    for (int j = 0; j < num_kmers; j++) {
        strncpy(kmer, read + j, kmer_len);

        exists = kmer_map.count(kmer) == 0 ? false : true;

        if (exists) {
            hash_key = kmer;
        } else {
            hash_key = new char [kmer_len + 1];
            strcpy(hash_key, kmer);
            kmer_keys.push_back(hash_key);
        }

        kmer_map[hash_key]++;
    }

    return retain;
}

inline int
kmer_lookup(SeqKmerHash &kmer_map,
            char *read, char *kmer,
            int num_kmers,
            int &rare_k, int &abundant_k)
{
    //
    // Generate kmers from this read, lookup kmer frequency in dataset.
    //
    rare_k         = 0;
    abundant_k     = 0;
    kmer[kmer_len] = '\0';
    int cnt        = 0;

    vector<int> cnts, sorted_cnts;
    cnts.reserve(num_kmers);
    sorted_cnts.reserve(num_kmers);

    // cout << "# " << read << "\n";

    for (int j = 0; j < num_kmers; j++) {
        strncpy(kmer, read + j, kmer_len);

        cnt = kmer_map[kmer];
        cnts.push_back(cnt);
        sorted_cnts.push_back(cnt);

        // cout << kmer << "\t" << j << "\t" << cnt << "\n";

        if (cnt >= max_k_freq) abundant_k++;
    }

    //
    // Calculate the median kmer frequency along the read.
    //
    sort(sorted_cnts.begin(), sorted_cnts.end());
    double median = num_kmers % 2 == 0 ?
        (double) (sorted_cnts[num_kmers / 2 - 1] + sorted_cnts[num_kmers / 2]) / 2.0 :
        (double) sorted_cnts[num_kmers / 2 - 1];
    // cout << "# median: " << median << "\n";

    double bound = round(median * min_k_pct);
    // cout << "# kmer cov bound: " << bound << "\n";

    //
    // Look for runs of rare kmers.
    //
    // We will slide a window across the read, f represents the front of the window, b
    // represents the back. Each  time a kmer is below the bound we will increment run_cnt,
    // which represents the number of kmers in the window below the bound. If 2/3 of the
    // kmers in the window go below the bound, assume a sequencing error has occurred.
    //
    int run_cnt = 0;
    int b = 0;
    for (int f = 0; f < num_kmers; f++) {
        if (f >= kmer_len) {
            b++;
            if (cnts[b] <= bound)
                run_cnt--;
        }
        if (cnts[f] <= bound) {
            run_cnt++;
            if (run_cnt >= min_lim) {
                rare_k++;
                // cout << "# Rejecting read, position: " << f << "; run_cnt: " << run_cnt << "\n";
                return 0;
            }
        }
        // cout << "# b: " << b << "; f: " << f << "; run_cnt: " << run_cnt << "; counts[front]: " << cnts[f] << "; bound: " << bound << "\n";
    }
    // cout << "\n\n";

    return 0;
}

// inline int
// kmer_lookup(SeqKmerHash &kmer_map,
//             char *read, char *kmer,
//             int num_kmers,
//             int &rare_k, int &abundant_k, bool &complex)
// {
//     //
//     // Generate kmers from this read, lookup kmer frequency in dataset.
//     //
//     rare_k         = 0;
//     abundant_k     = 0;
//     complex        = false;
//     kmer[kmer_len] = '\0';
//     int cnt        = 0;
//     int rare_k_lim = (int) round((double) kmer_len * (1.0/2.0));

//     vector<int> cnts;
//     cnts.reserve(num_kmers);

//     // cout << "# " << read << "\n";

//     for (int j = 0; j < num_kmers; j++) {
//         strncpy(kmer, read + j, kmer_len);

//         cnt = kmer_map[kmer];

//         if (cnt >= 100000)
//             cnts.push_back(100000);
//         else if (cnt >= 10000)
//             cnts.push_back(10000);
//         else if (cnt >= 1000)
//             cnts.push_back(1000);
//         else if (cnt >= 100)
//             cnts.push_back(100);
//         else if (cnt >= 10)
//             cnts.push_back(10);
//         else
//             cnts.push_back(1);

//         // cout << kmer << "\t" << j << "\t" << cnt << "\n";

//         if (cnt >= max_k_freq) abundant_k++;
//     }

//     // //
//     // // Detect the number of kmer coverage transitions.
//     // //
//     // int t   = 0;
//     // int cov = cnts[0];
//     // cout << "\nDetermining transitions:\n" << kmer << "\t" << "0" << "\t" << cnts[0] << "\n";
//     // for (int j = 1; j < num_kmers; j++)
//     //         if (cnts[j] != cov) {
//     //             cov = cnts[j];
//     //             t++;
//     //             cout << kmer << "\t" << j << "\t" << cnts[j] << ": Transition." << "\n";
//     //         } else {
//     //             cout << kmer << "\t" << j << "\t" << cnts[j] << "\n";
//     //         }
//     // cerr << t << " total cnts.\n";

//     //
//     // Look for runs of kmers at various orders of magnitude.
//     //
//     // We will slide a window across the read, f represents the front of the window, b
//     // represents the back. Each  time a kmer is below the bound we will increment run_cnt,
//     // which represents the number of kmers in the window below the bound. If 2/3 of the
//     // kmers in the window go below the bound, assume a sequencing error has occurred.

//     // Run counters:
//     //   1      10      100      1k     10k    100k
//     // runs[0] runs[1] runs[2] runs[3] runs[4] runs[5]
//     int runs[6] = {0};
//     int prev_cnt, run_cnt, tot_trans;
//     int f = 0;

//     while (f < num_kmers) {

//         tot_trans = 0;
//         run_cnt   = 1;
//         prev_cnt  = cnts[f];
//         f++;

//         while (f < num_kmers && cnts[f] == prev_cnt) {
//             // cout << "# window front: " << f << "; run_cnt: " << run_cnt << "; prev_cnt: " << prev_cnt << "\n";
//             f++;
//             run_cnt++;
//         }

//         if (run_cnt >= rare_k_lim) {
//             // cout << "#   found transition run, position: " << f-1 << "; run_cnt: " << run_cnt << "\n";
//             switch(prev_cnt) {
//             case 1:
//                 runs[0]++;
//                 break;
//             case 10:
//                 runs[1]++;
//                 break;
//             case 100:
//                 runs[2]++;
//                 break;
//             case 1000:
//                 runs[3]++;
//                 break;
//             case 10000:
//                 runs[4]++;
//                 break;
//             case 100000:
//                 runs[5]++;
//                 break;
//             }
//         }

//         for (int j = 0; j < 6; j++)
//             if (runs[j] > 0) tot_trans++;

//         // cout << "# Total transitions: " << tot_trans << "\n";

//         if (tot_trans >= transition_lim) {
//             // cout << "# Rejecting read.\n";
//             rare_k++;
//             return 0;
//         }
//     }

//     // cout << "\n\n";

//     return 0;
// }

int
free_kmer_hash(SeqKmerHash &kmer_map, vector<char *> &kmer_map_keys)
{
    for (uint i = 0; i < kmer_map_keys.size(); i++) {
        delete [] kmer_map_keys[i];
    }
    kmer_map_keys.clear();

    kmer_map.clear();

    return 0;
}

int print_results(map<string, map<string, long> > &counters) {
    map<string, map<string, long> >::iterator it;

    string log_path = out_path + "kmer_filter.log";
    ofstream log(log_path.c_str());

    if (log.fail()) {
        cerr << "Unable to open log file '" << log_path << "'\n";
        return 0;
    }

    cerr << "Outputing details to log: '" << log_path << "'\n\n";

    log << "File\t"
        << "Retained Reads\t"
        << "Rare K\t"
        << "Abundant K\t"
        << "Total\n";

    for (it = counters.begin(); it != counters.end(); it++) {
        log << it->first                 << "\t"
            << it->second["retained"]    << "\t"
            << it->second["rare_k"]      << "\t"
            << it->second["abundant_k"]  << "\t"
            << it->second["total"]       << "\n";
    }

    map<string, long> c;
    c["total"]       = 0;

    //
    // Total up the individual counters
    //
    for (it = counters.begin(); it != counters.end(); it++) {
        c["total"]       += it->second["total"];
        c["retained"]    += it->second["retained"];
        c["rare_k"]      += it->second["rare_k"];
        c["abundant_k"]  += it->second["abundant_k"];
    }

    cerr <<
        c["total"] << " total sequences;\n"
         << "  " << c["rare_k"]      << " rare k-mer reads;\n"
         << "  " << c["abundant_k"]  << " abundant k-mer reads;\n"
         << c["retained"] << " retained reads.\n";

    log        << "Total Sequences\t"      << c["total"]       << "\n"
        << "Retained Reads\t"       << c["retained"]      << "\n";

    log.close();

    return 0;
}

int build_file_list(vector<string> &in_files, vector<pair<string, string> > &files) {
    string file, suffix;
    int    pos;

    //
    // Scan a directory for a list of files.
    //
    if (in_path.length() > 0) {
        struct dirent *direntry;

        DIR *dir = opendir(in_path.c_str());

        if (dir == NULL) {
            cerr << "Unable to open directory '" << in_path << "' for reading.\n";
            exit(1);
        }

        while ((direntry = readdir(dir)) != NULL) {
            file = direntry->d_name;

            if (file.substr(0, 1) == ".")
                continue;

            //
            // If the file is gzip'ed, remove the '.gz' suffix.
            //
            pos  = file.find_last_of(".");
            if ((in_file_type == FileT::gzfastq || in_file_type == FileT::gzfasta) &&
                file.substr(pos) == ".gz") {
                file = file.substr(0, pos);
                pos  = file.find_last_of(".");
            }

            //
            // Check that the remaining file name has the right suffix.
            //
            suffix = file.substr(pos + 1);
            if (in_file_type == FileT::fastq && (suffix.substr(0, 2) == "fq" || suffix.substr(0, 5) == "fastq"))
                files.push_back(make_pair(in_path, file));
            else if (in_file_type == FileT::fasta && (suffix.substr(0, 2) == "fa" || suffix.substr(0, 5) == "fasta"))
                files.push_back(make_pair(in_path, file));
        }

        if (files.size() == 0)
            cerr << "Unable to locate any input files to process within '" << in_path << "'\n";

    } else {
        string path;

        for (uint i = 0; i < in_files.size(); i++) {
            //
            // Files specified directly:
            //    Break off file path and store path and file name.
            //    Check if this is a gzip'ed file and if so, remove 'gz' suffix.
            //
            file = in_files[i];
            pos  = file.find_last_of(".");
            if ((in_file_type == FileT::gzfastq || in_file_type == FileT::gzfasta) &&
                file.substr(pos) == ".gz") {
                file = file.substr(0, pos);
                pos  = file.find_last_of(".");
            }
            pos  = file.find_last_of("/");
            path = file.substr(0, pos + 1);
            files.push_back(make_pair(path, file.substr(pos+1)));
        }
    }

    return 0;
}

int parse_command_line(int argc, char* argv[]) {
    string pair_1, pair_2;
    int c;

    while (1) {
        static struct option long_options[] = {
            {"help",         no_argument,       NULL, 'h'},
            {"version",      no_argument,       NULL, 'v'},
            {"discards",     no_argument,       NULL, 'D'},
            {"pair-1",       required_argument, NULL, '1'}, {"pair_1",       required_argument, NULL, '1'},
            {"pair-2",       required_argument, NULL, '2'}, {"pair_2",       required_argument, NULL, '2'},
            {"infile-type",  required_argument, NULL, 'i'}, {"infile_type",  required_argument, NULL, 'i'},
            {"outfile-type", required_argument, NULL, 'y'}, {"outfile_type", required_argument, NULL, 'y'},
            {"file",         required_argument, NULL, 'f'},
            {"path",         required_argument, NULL, 'p'},
            {"outpath",      required_argument, NULL, 'o'},
            {"k-dist",       no_argument,       NULL, 'I'}, {"k_dist",       no_argument,       NULL, 'I'},
            {"rare",         no_argument,       NULL, 'R'},
            {"abundant",     no_argument,       NULL, 'A'},
            {"normalize",    required_argument, NULL, 'N'},
            {"k-len",        required_argument, NULL, 'K'}, {"k_len",        required_argument, NULL, 'K'},
            {"max-k-freq",   required_argument, NULL, 'M'}, {"max_k_freq",   required_argument, NULL, 'M'},
            {"min-lim",      required_argument, NULL, 'F'}, {"min_lim",      required_argument, NULL, 'F'},
            {"max-lim",      required_argument, NULL, 'G'}, {"max_lim",      required_argument, NULL, 'G'},
            {"min-k-pct",    required_argument, NULL, 'P'}, {"min_k_pct",    required_argument, NULL, 'P'},
            {"read-k-freq",  required_argument, NULL, 'r'}, {"read_k_freq",  required_argument, NULL, 'r'},
            {"write-k-freq", required_argument, NULL, 'w'}, {"write_k_freq", required_argument, NULL, 'w'},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;

        c = getopt_long(argc, argv, "hvRADkP:N:I:w:r:K:F:G:M:m:i:y:f:o:t:p:1:2:", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            help();
            break;
        case 'i':
             if (strcasecmp(optarg, "fasta") == 0)
                in_file_type = FileT::fasta;
             else if (strcasecmp(optarg, "gzfasta") == 0)
                 in_file_type = FileT::gzfasta;
             else if (strcasecmp(optarg, "gzfastq") == 0)
                 in_file_type = FileT::gzfastq;
             else
                 in_file_type = FileT::fastq;
            break;
        case 'y':
            if (strcasecmp(optarg, "fasta") == 0)
                out_file_type = FileT::fasta;
            else
                out_file_type = FileT::fastq;
            break;
        case 'f':
            in_files.push_back(optarg);
            break;
        case '1':
            pair_1 = optarg;
            break;
        case '2':
            pair_2 = optarg;
            if (pair_1.length() == 0) help();
            in_pair_files.push_back(pair_1);
            in_pair_files.push_back(pair_2);
            pair_1 = "";
            pair_2 = "";
            break;
        case 'p':
            in_path = optarg;
            break;
        case 'o':
            out_path = optarg;
            break;
        case 'D':
            discards = true;
            break;
        case 'I':
            kmer_distr = true;
            break;
        case 'R':
            filter_mare_k = true;
            break;
        case 'A':
            filter_abundant_k = true;
            break;
        case 'N':
            normalize = true;
            normalize_lim = is_integer(optarg);
            break;
        case 'K':
            kmer_len = is_integer(optarg);
            break;
        case 'M':
            max_k_freq = is_integer(optarg);
            break;
        case 'F':
            min_lim = is_integer(optarg);
            break;
        case 'G':
            max_lim = is_integer(optarg);
            break;
        case 'P':
            min_k_pct = is_double(optarg);
            break;
        case 'r':
            read_k_freq = true;
            k_freq_path = optarg;
            break;
        case 'w':
            write_k_freq = true;
            k_freq_path = optarg;
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

    if (in_files.size() == 0 && in_pair_files.size() == 0 && in_path.length() == 0) {
        cerr << "You must specify an input file of a directory path to a set of input files.\n";
        help();
    }

    if (in_files.size() > 0 && in_path.length() > 0) {
        cerr << "You must specify either a single input file (-f) or a directory path (-p), not both.\n";
        help();
    }

    if (in_path.length() > 0 && in_path.at(in_path.length() - 1) != '/')
        in_path += "/";

    if (out_path.length() == 0)
        out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/')
        out_path += "/";

    if (in_file_type == FileT::unknown)
        in_file_type = FileT::fastq;

    if (read_k_freq && write_k_freq) {
        cerr << "You may either read a set of kmer frequencies, or write kmer frequencies, not both.\n";
        help();
    }

    if (min_k_pct < 0.0 || min_k_pct > 1.0) {
        cerr << "Percentage to consider a kmer rare must be between 0 and 1.0.\n";
        help();
    }

    //
    // Check that the output path exists.
    //
    struct stat info;
    if (stat(out_path.c_str(), &info) != 0) {
        cerr << "Unable to locate the specified output path, '" << out_path << "'\n";
        exit(1);
    }

    return 0;
}

void version() {
    cerr << "kmer_filter " << VERSION << "\n\n";

    exit(1);
}

void help() {
    cerr << "kmer_filter " << VERSION << "\n"
              << "kmer_filter [-f in_file_1 [-f in_file_2...] | -p in_dir] [-1 pair_1 -2 pair_2 [-1 pair_1...]] -o out_dir [-i type] [-y type] [-D] [-h]\n"
              << "  f: path to the input file if processing single-end seqeunces.\n"
              << "  i: input file type, either 'bustard' for the Illumina BUSTARD output files, 'fasta', 'fastq', 'gzfasta', or 'gzfastq' (default 'fastq').\n"
              << "  p: path to a directory of files (for single-end files only).\n"
              << "  1: specify the first in a pair of files to be processed together.\n"
              << "  2: specify the second in a pair of files to be processed together.\n"
              << "  o: path to output the processed files.\n"
              << "  y: output type, either 'fastq' or 'fasta' (default fastq).\n"
              << "  D: capture discarded reads to a file.\n"
              << "  h: display this help messsage.\n\n"
              << "  Filtering options:\n"
              << "    --rare: turn on filtering based on rare k-mers.\n"
              << "    --abundant: turn on filtering based on abundant k-mers.\n"
              << "    --k-len <len>: specify k-mer size (default 15).\n\n"
              << "  Advanced filtering options:\n"
              << "    --max-k-freq <value>: specify the number of times a kmer must occur to be considered abundant (default 20,000).\n"
              << "    --min-lim <value>: specify number of rare kmers occuring in a row required to discard a read (default 80% of the k-mer length).\n"
              << "    --max-lim <value>: specify number of abundant kmers required to discard a read (default 80% of the k-mers in a read).\n\n"
              << "  Normalize data:\n"
              << "    --normalize <depth>: normalize read depth according to k-mer coverage.\n\n"
              << "  Characterizing K-mers:\n"
              << "    --write-k-freq: write kmers along with their frequency of occurrence and exit.\n"
              << "    --k-dist: print k-mer frequency distribution and exit.\n\n"
              << "  Advanced input options:\n"
              << "    --read-k-freq <path>: read a set of kmers along with their frequencies of occurrence instead of reading raw input files.\n"
              << "\n";

    exit(1);
}
