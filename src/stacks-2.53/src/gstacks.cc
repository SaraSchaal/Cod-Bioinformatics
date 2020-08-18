// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2017-2018, Julian Catchen <jcatchen@illinois.edu>
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
#include <getopt.h>
#include <zlib.h>

#include "gstacks.h"

#include "htslib/bgzf.h"
#include "constants.h"
#include "utils.h"
#include "log_utils.h"
#include "catalog_utils.h"
#include "locus.h"
#include "locus_readers.h"
#include "debruijn.h"
#include "aln_utils.h"
#include "Alignment.h"
#include "models.h"

//
// Argument globals.
//
GStacksInputT  input_type = GStacksInputT::unknown;
string         popmap_path;  // Set if --popmap is given. For logging purposes.
vector<string> sample_names; // Set if --popmap is given.
vector<string> in_bams;
vector<string> out_bams;

string out_dir;

BamCLocBuilder::Config refbased_cfg {true, 1000, 10, 0.20, 1, false, 1};

int    num_threads        = 1;
bool   quiet              = false;
bool   ignore_pe_reads    = false;
double min_aln_cov        = 0.75;
int    min_se_pe_overlap  = 5;
double overlap_min_pct_id = 0.80;
bool   bam_output         = false;
bool   detailed_output    = false;
bool   rm_unpaired_reads  = false;
bool   rm_pcr_duplicates  = false;

modelt model_type = marukilow;
unique_ptr<const Model> model;

size_t km_length         = 31;
size_t max_debruijn_reads = 1000;
size_t min_km_count      = 2;
size_t max_fragment_alns = 2;

pair<size_t,size_t> phasing_cooccurrences_thr_range = {2, 2};
bool phasing_dont_prune_hets = false;
long phasing_min_mac = 1;

bool   dbg_no_overlaps     = false;
bool   dbg_no_haplotypes   = false;
bool   dbg_print_cloc_ids  = false;
bool   dbg_write_gfa       = false;
bool   dbg_write_gfa_notdag = false;
bool   dbg_write_alns      = false;
bool   dbg_write_hapgraphs = false;
bool   dbg_write_nt_depths = false;
bool   dbg_log_stats_phasing = false;
bool   dbg_phasing_no_2ndpass = false;
size_t dbg_denovo_min_loc_samples = 0;

//
// Additional globals.
//
const string prog_name = "gstacks";
unique_ptr<LogAlterator> logger;
BGZF* o_gzfasta_f = NULL;
unique_ptr<VcfWriter> o_vcf_f;
unique_ptr<VersatileWriter> o_details_f;

ofstream o_aln_f;
ofstream o_hapgraphs_f;
const char o_aln_header[] =
    "# This prints observed read haplotypes:\n"
    "# show_loc() { loc=$1; cat ./gstacks.alns | sed -n \"/^END $loc\\b/ q; /^BEGIN $loc\\b/,$ p\" | tail -n+2; }\n"
    "# snp_cols() { loc=$1; zcat ./gstacks.calls | awk \"\\$1==$loc; \\$1>$loc {exit}\" | awk '$5!=\".\"' | cut -f2 | paste -sd ','; }\n"
    "# show_haps() { loc=$1; cols=$2; spl=$3; show_loc $loc | grep \"\\b$spl\\b\" | cut -f3 | cut -c \"${cols//./,}\" | sort; }\n"
    "# format_haps() { loc=$1; cols=$2; spl=$3; show_haps $loc $cols $spl | transpose-ch | paste <(echo $cols | tr ,. '\\n' | xargs printf '% 4s\\n') - | sed 's/\\t/  /' | seq.color | { echo; echo $loc/$spl; echo; cat; echo; }; }\n"
    "# true_loci() { loc=$1; spl=$2; show_loc $loc | grep \"\\b$spl\\b\" | grep -v ref | cut -d: -f1 | sort -u; }\n"
    ;
const char o_hapgraphs_header[] =
    "# dot -Tpdf -O gstacks.hapgraphs.dot\n"
    "# loc=371\n"
    "# { g=gstacks.hapgraphs.dot; sed -n '0,/^subgraph/p' $g | head -n-1; sed -n \"/^subgraph cluster_loc$loc\\b/,/^}/p\" $g; echo \\}; } | dot -Tpdf -o haps.$loc.pdf\n"
    "graph {\n"
    "edge[color=\"grey60\",fontsize=12,labeljust=\"l\"];\n";
    ;

//
// main
// ==========
//

int
main(int argc, char** argv)
{
try {
    //
    // Parse arguments.
    //
    string opts_report = parse_command_line(argc, argv);

    //
    // Open the BAM file(s).
    //
    vector<Bam*> bam_f_ptrs;
try {
    for (const string& in_bam : in_bams)
        bam_f_ptrs.push_back(new Bam(in_bam));
} catch(exception& e) {
    if (bam_f_ptrs.size() >= 250)
        cerr << "Error: You might need to increase your system's max-open-files limit,"
                " see https://groups.google.com/d/msg/stacks-users/GZqJM_WkMhI/m9Hreg4oBgAJ\n";
    throw;
}

    //
    // Open the log.
    //
    logger.reset(new LogAlterator(out_dir + "gstacks", true, quiet, argc, argv));
    cout << opts_report << flush;

    //
    // Initialize the locus readers.
    //
    cout << "\nReading BAM headers...\n" << flush;
    unique_ptr<BamCLocReader> bam_cloc_reader;
    unique_ptr<BamCLocBuilder> bam_cloc_builder;
    const MetaPopInfo* bam_mpopi;
    if (input_type == GStacksInputT::denovo_popmap || input_type == GStacksInputT::denovo_merger) {
        bam_cloc_reader.reset(new BamCLocReader(move(bam_f_ptrs)));
        // (Note: For de novo (Stacks) input we use read groups, even when a
        // popmap was specified, since the input files ought to have them and
        // this allows to detect corruption.)
        bam_mpopi = &bam_cloc_reader->mpopi();
    } else if (input_type == GStacksInputT::refbased_popmap || input_type == GStacksInputT::refbased_list) {
        bam_cloc_builder.reset(new BamCLocBuilder(move(bam_f_ptrs), refbased_cfg, sample_names));
        bam_mpopi = &bam_cloc_builder->mpopi();
    } else {
        DOES_NOT_HAPPEN;
    }

    //
    // Open the output files.
    //
    string o_gzfasta_path = out_dir + "catalog.fa.gz";
    o_gzfasta_f = bgzf_open(o_gzfasta_path.c_str(), "w");
    check_open(o_gzfasta_f, o_gzfasta_path);

    VcfHeader vcf_header;
    vcf_header.add_std_meta();
    for(size_t s : bam_mpopi->sample_indexes_orig_order())
        vcf_header.add_sample(bam_mpopi->samples()[s].name);
    o_vcf_f.reset(new VcfWriter(out_dir + "catalog.calls", move(vcf_header)));

    if (detailed_output) {
        o_details_f.reset(new VersatileWriter(out_dir + "gstacks.details.gz"));
        *o_details_f << "# show_loc() { loc=$1; zcat ./gstacks.details.gz | sed -rn \"/^BEGIN locus $loc\\b/,\\$ p; /^END locus $loc\\b/ q;\"; }\n";
    }

    vector<unique_ptr<Bam>> bam_of_ptrs;
    if (bam_output) {
    try {
        stringstream header_ss;
        header_ss << "@HD\tVN:1.5\tSO:coordinate\n";
        for (const Sample& s : bam_cloc_reader->mpopi().samples())
            header_ss << "@RG\tID:" << s.id << "\tSM:" << s.name << "\tid:" << s.id << '\n';
        for (size_t target_i=0; target_i<bam_cloc_reader->n_loci(); ++target_i)
            header_ss << "@SQ\tSN:" << bam_cloc_reader->target2id(target_i) << "\tLN:10000\n";
        BamHeader header = BamHeader(header_ss.str());
        if (input_type == GStacksInputT::denovo_merger) {
            assert(out_bams.size() == 1);
            bam_of_ptrs.emplace_back(new Bam(out_bams[0], BamHeader(header)));
        } else if (input_type == GStacksInputT::denovo_popmap) {
            assert(out_bams.size() == bam_cloc_reader->bam_fs().size());
            for(size_t i=0; i<out_bams.size(); ++i) {
                bam_of_ptrs.emplace_back(new Bam(out_bams[i], BamHeader(header)));
            }
        } else {
            DOES_NOT_HAPPEN;
        }
    } catch(exception& e) {
        if (bam_f_ptrs.size() + bam_of_ptrs.size() >= 250)
            cerr << "Error: You might need to increase your system's max-open-files limit,"
                    " see https://groups.google.com/d/msg/stacks-users/GZqJM_WkMhI/m9Hreg4oBgAJ\n";
        throw;
    }}

    if (dbg_write_alns) {
        string o_aln_path = out_dir + "gstacks.alns";
        o_aln_f.open(o_aln_path);
        check_open(o_aln_f, o_aln_path);
        o_aln_f << o_aln_header;
    }

    if (dbg_write_hapgraphs) {
        string o_hapgraphs_path = out_dir + "gstacks.hapgraphs.dot";
        o_hapgraphs_f.open(o_hapgraphs_path);
        check_open(o_hapgraphs_f, o_hapgraphs_path);
        o_hapgraphs_f << o_hapgraphs_header;
    }

    //
    // Process every locus
    //
    GenotypeStats gt_stats (bam_mpopi->samples().size());
    HaplotypeStats hap_stats (bam_mpopi->samples().size());
    ContigStats denovo_ctg_stats {};

    // For clocking.
    Timers t_threads_totals;
    Timer t_parallel;
    Timer t_writing_vcf;
    size_t n_writes = 0;
    size_t max_size_before_write = 0;

    // For parallelization.
    int omp_return = 0;
    std::deque<pair<bool,string>> fa_outputs;
    std::deque<pair<bool,string>> vcf_outputs;
    std::deque<pair<bool,string>> det_outputs; // (details)
    std::deque<pair<bool,vector<vector<BamRecord>>>> bam_outputs;
    size_t next_fa_to_write = 0; // locus index, in the input BAM file.
    size_t next_vcf_to_write = 0;
    size_t next_det_to_write = 0;
    size_t next_bam_to_write = 0;

    cout << "Processing all loci...\n" << flush;
    ProgressMeter progress =
            (input_type == GStacksInputT::denovo_popmap || input_type == GStacksInputT::denovo_merger) ?
            ProgressMeter(cout, true, bam_cloc_reader->tally_n_components())
            : ProgressMeter(cout, false, 1000);

    bool eof = false;
    #pragma omp parallel num_threads(num_threads)
    { try {
        LocusProcessor loc_proc (bam_mpopi->samples().size(), bam_cloc_reader.get());
        Timers& t = loc_proc.timers();

        CLocReadSet loc (*bam_mpopi); // For denovo.
        CLocAlnSet aln_loc; // For ref-based.
        bool thread_eof = false;
        while(omp_return == 0) {
            t.reading.restart();
            #pragma omp critical(read)
            { try {
                if (!eof && omp_return == 0) {
                    if (input_type == GStacksInputT::denovo_popmap || input_type == GStacksInputT::denovo_merger)
                        eof = !bam_cloc_reader->read_one_locus(loc);
                    else
                        eof = !bam_cloc_builder->build_one_locus(aln_loc);
                }
                thread_eof = eof;
            } catch (exception& e) {
                omp_return = stacks_handle_exceptions(e);
            }}
            t.reading.update();
            if (thread_eof || omp_return != 0)
                break;

            // Process it.
            t.processing.restart();
            size_t loc_i;
            if (input_type == GStacksInputT::denovo_popmap || input_type == GStacksInputT::denovo_merger) {
                loc_i = loc.bam_i();
                if (dbg_print_cloc_ids)
                    cerr << (to_string(loc.id()) + '\n') << flush;
                if (ignore_pe_reads)
                    loc.pe_reads().clear();
                loc_proc.process(loc);
            } else {
                assert(aln_loc.id() >= 1);
                loc_i = aln_loc.id() - 1;
                if (dbg_print_cloc_ids)
                    cerr << (to_string(aln_loc.id()) + '\n') << flush;
                if (refbased_cfg.paired)
                    aln_loc.merge_paired_reads();
                loc_proc.process(aln_loc);
            }
            t.processing.update();

            // Write the FASTA output.
            t.writing_fa.restart();
            #pragma omp critical(write_fa)
            {
                for (size_t i=next_fa_to_write+fa_outputs.size(); i<=loc_i; ++i)
                    fa_outputs.push_back( {false, string()} );
                fa_outputs[loc_i - next_fa_to_write] = {true, move(loc_proc.fasta_out())};

                while (!fa_outputs.empty() && fa_outputs.front().first) {
                    const string& fa = fa_outputs.front().second;
                    if (!fa.empty() && bgzf_write(o_gzfasta_f, fa.c_str(), fa.length()) <= 0)
                        throw std::ios::failure("gzwrite");
                    fa_outputs.pop_front();
                    ++next_fa_to_write;
                }
            }
            t.writing_fa.update();

            // Write the VCF output.
            t.writing_vcf.restart();
            #pragma omp critical(write_vcf)
            {
                for (size_t i=next_vcf_to_write+vcf_outputs.size(); i<=loc_i; ++i)
                    vcf_outputs.push_back( {false, string()} );
                vcf_outputs[loc_i - next_vcf_to_write] = {true, move(loc_proc.vcf_out()) };

                if (vcf_outputs.front().first) {
                    ++n_writes;
                    if (vcf_outputs.size() > max_size_before_write)
                        max_size_before_write = vcf_outputs.size();
                    t_writing_vcf.restart();
                    do {
                        o_vcf_f->file() << vcf_outputs.front().second;
                        vcf_outputs.pop_front();
                        if (input_type == GStacksInputT::denovo_popmap || input_type == GStacksInputT::denovo_merger)
                            progress += bam_cloc_reader->n_catalog_components_of(next_vcf_to_write);
                        else
                            ++progress;
                        ++next_vcf_to_write;
                    } while (!vcf_outputs.empty() && vcf_outputs.front().first);
                    t_writing_vcf.update();
                }
            }
            t.writing_vcf.update();

            // Write the alignments.
            if (bam_output) {
                t.writing_bams.restart();
                #pragma omp critical(write_bams)
                {
                    for (size_t i=next_bam_to_write+bam_outputs.size(); i<=loc_i; ++i)
                        bam_outputs.push_back( make_pair(false, vector<vector<BamRecord>>()) );
                    bam_outputs[loc_i - next_bam_to_write] = make_pair(true, move(loc_proc.bam_out()));

                    while (!bam_outputs.empty() && bam_outputs.front().first) {
                        vector<vector<BamRecord>>& recs = bam_outputs.front().second;
                        if (!recs.empty()) {
                            assert(recs.size() == bam_of_ptrs.size());
                            for (size_t i=0; i<recs.size(); ++i)
                                for (const BamRecord& r : recs[i])
                                    bam_of_ptrs[i]->write(r);
                        }
                        bam_outputs.pop_front();
                        ++next_bam_to_write;
                    }
                }
                t.writing_bams.update();
            }

            // Write the detailed output.
            if (detailed_output) {
                t.writing_details.restart();
                #pragma omp critical(write_details)
                {
                    for (size_t i=next_det_to_write+det_outputs.size(); i<=loc_i; ++i)
                        det_outputs.push_back( {false, string()} );
                    det_outputs[loc_i - next_det_to_write] = {true, move(loc_proc.details_out())};

                    while (!det_outputs.empty() && det_outputs.front().first) {
                        *o_details_f << det_outputs.front().second;
                        det_outputs.pop_front();
                        ++next_det_to_write;
                    }
                }
                t.writing_details.update();
            }
        }

        // Tally the per-thread statistics.
        #pragma omp critical(stats)
        if (omp_return == 0) {
            denovo_ctg_stats += loc_proc.ctg_stats();
            gt_stats  += loc_proc.gt_stats();
            hap_stats  += loc_proc.hap_stats();
            t_threads_totals += loc_proc.timers();
        }
    } catch (exception& e) {
        #pragma omp critical(exc)
        omp_return = stacks_handle_exceptions(e);
    }}
    if (omp_return != 0)
        return omp_return;
    t_parallel.update();
    progress.done();

    //
    // Report statistics.
    //
    ostream o_fp1 (cout.rdbuf());
    o_fp1 << std::fixed << std::setprecision(1);
    ostream o_fp3 (cout.rdbuf());
    o_fp3 << std::fixed << std::setprecision(3);
    ostream x_fp1 (logger->x.rdbuf());
    x_fp1 << std::fixed << std::setprecision(1);
    ostream x_fp3 (logger->x.rdbuf());
    x_fp3 << std::fixed << std::setprecision(3);

    if (input_type == GStacksInputT::denovo_popmap || input_type == GStacksInputT::denovo_merger) {
        if (denovo_ctg_stats.n_loci_w_pe_reads == 0) {
            cout << "Input appears to be single-end (no paired-end reads were seen).\n";
        } else {
            // Report assembly statistics.
            assert(!ignore_pe_reads);
            const ContigStats& cs = denovo_ctg_stats;
            size_t no_pe   = cs.n_loci_no_pe_reads() + cs.n_loci_almost_no_pe_reads;
            size_t pe_ndag  = cs.n_loci_pe_graph_not_dag;
            size_t pe_ctg = cs.n_loci_ctg();
            auto pct = [&cs](size_t n) { return as_percentage((double) n / cs.n_nonempty_loci); };

            o_fp1 << "\n"
               << "Attempted to assemble and align paired-end reads for " << cs.n_nonempty_loci << " loci:\n"
               << "  " << no_pe << " loci had no or almost no paired-end reads (" << pct(no_pe) << ");\n"
               << "  " << pe_ndag << " loci had paired-end reads that couldn't be assembled into a contig ("
               << pct(pe_ndag) << ");\n"
               << "  For the remaining " << pe_ctg << " loci (" << pct(pe_ctg) << "), a paired-end contig was assembled;\n"
               << "    Average contig size was " << cs.ctg_avg_length() << " bp;\n"
               << "  " << cs.n_overlaps << " paired-end contigs overlapped the forward region ("
               << as_percentage((double) cs.n_overlaps / cs.n_loci_ctg()) << ")\n"
               << "    Mean overlap: " << cs.mean_olap_length() << "bp; mean size of overlapped loci after merging: "
               << cs.mean_olapd_locus_length() << ";\n"
               << "  Out of " << cs.n_tot_reads << " paired-end reads in these loci (mean "
               << (double) cs.n_aln_reads / pe_ctg << " reads per locus),\n"
               << "    " << cs.n_aln_reads << " were successfuly aligned ("
               << as_percentage((double) cs.n_aln_reads / cs.n_tot_reads) << ");\n"
               << "  Mean insert length was " << cs.insert_length_olap_mv.mean() << ", stdev: "
               << cs.insert_length_olap_mv.sd_p() << " (based on aligned reads in overlapped loci).\n";
        }

    } else if (input_type == GStacksInputT::refbased_popmap || input_type == GStacksInputT::refbased_list) {
        // Report statistics on the input BAM(s).
        BamCLocBuilder::BamStats bam_stats {};
        const vector<BamCLocBuilder::BamStats>& bam_stats_s = bam_cloc_builder->bam_stats_per_sample();
        for (auto& bstats : bam_stats_s)
            bam_stats += bstats;
        size_t tot = bam_stats.n_primary + bam_stats.n_unmapped;
        cout << "\n"
             << "Read " << bam_stats.n_records << " BAM records:\n"
             << "  kept " << bam_stats.n_primary_kept() << " primary alignments ("
             << as_percentage(bam_stats.n_primary_kept(), tot) << ")";
        if (refbased_cfg.paired)
            cout << ", of which " << bam_stats.n_primary_kept_read2 << " reverse reads";
        cout << "\n"
             << "  skipped " << bam_stats.n_primary_mapq << " primary alignments with insufficient mapping qualities ("
             << as_percentage((double) bam_stats.n_primary_mapq / tot) << ")\n"
             << "  skipped " << bam_stats.n_primary_softclipped << " excessively soft-clipped primary alignments ("
             << as_percentage((double) bam_stats.n_primary_softclipped / tot) << ")\n"
             << "  skipped " << bam_stats.n_unmapped << " unmapped reads ("
             << as_percentage((double) bam_stats.n_unmapped / tot) << ")\n";
        if (bam_stats.n_secondary > 0 || bam_stats.n_supplementary > 0)
            cout << "  skipped some suboptimal (secondary/supplementary) alignment records\n";
        if (refbased_cfg.ign_pe_reads)
            cout << "  ignored " << bam_stats.n_ignored_read2_recs << " READ2 records\n";

        // Per sample BAM stats.
        logger->x << "\n"
                  << "BEGIN bam_stats_per_sample\n"
                  << "sample"
                     "\trecords"
                     "\tprimary_kept"
                     "\tkept_frac"
                     "\tprimary_kept_read2"
                     "\tprimary_disc_mapq"
                     "\tprimary_disc_sclip"
                     "\tunmapped"
                     "\tsecondary"
                     "\tsupplementary\n";
        size_t min_recs = SIZE_MAX;
        size_t max_recs = 0;
        double mean_recs = 0.0;
        double min_keep_rate = DBL_MAX;
        double max_keep_rate = 0.0;
        assert(bam_stats_s.size() == bam_mpopi->samples().size());
        for (size_t sample=0; sample<bam_stats_s.size(); ++sample) {
            auto& bstats = bam_stats_s[sample];
            x_fp3 << bam_mpopi->samples()[sample].name << '\t'
                  << bstats.n_records << '\t'
                  << bstats.n_primary_kept() << '\t'
                  << (double) bstats.n_primary_kept() / bstats.n_records << '\t'
                  << bstats.n_primary_kept_read2 << '\t'
                  << bstats.n_primary_mapq << '\t'
                  << bstats.n_primary_softclipped << '\t'
                  << bstats.n_unmapped << '\t'
                  << bstats.n_secondary << '\t'
                  << bstats.n_supplementary << '\n';
            mean_recs += bstats.n_records;
            if (bstats.n_records < min_recs)
                min_recs = bstats.n_records;
            if (bstats.n_records > max_recs)
                max_recs = bstats.n_records;
            double keep_rate = (double) bstats.n_primary_kept() / bstats.n_records;
            if (keep_rate < min_keep_rate)
                min_keep_rate = keep_rate;
            if (keep_rate > max_keep_rate)
                max_keep_rate = keep_rate;
        }
        logger->x << "END bam_stats_per_sample\n";
        mean_recs /= bam_mpopi->samples().size();
        o_fp1 << "\n"
              << "  Per-sample stats (details in '"
              << logger->distribs_path.substr(logger->distribs_path.rfind('/') + 1)
              << "'):\n"
              << "    read " << mean_recs << " records/sample (" << min_recs
              << "-" << max_recs << ")\n"
              << "    kept " << as_percentage(min_keep_rate) << "-"
              << as_percentage(max_keep_rate) << " of these\n";

        // Report statistics on the loci that were built.
        const BamCLocBuilder::LocStats& loc_stats = bam_cloc_builder->loc_stats();
        cout << "\n"
             << "Built " << loc_stats.n_loci_built << " loci comprising "
             << loc_stats.n_fw_reads;
        if (refbased_cfg.paired) {
            o_fp1 << " forward reads and "
                  << loc_stats.n_read_pairs() << " matching paired-end reads;"
                  << " mean insert length was " << loc_stats.insert_lengths_mv.mean()
                  << " (sd: " << loc_stats.insert_lengths_mv.sd_p() << ").\n";
        } else {
            cout << " reads.\n";
        }
    }

    if (rm_unpaired_reads) {
        // Report statistics on pairs.
        size_t n_unpaired = gt_stats.n_unpaired_reads_rm();
        size_t n_pcr_dupl = gt_stats.n_read_pairs_pcr_dupl();
        size_t n_used = gt_stats.n_read_pairs_used();
        size_t n_remaining_loci = gt_stats.n_genotyped_loci;
        cout << "Removed " << n_unpaired << " unpaired (forward) reads ("
             << as_percentage((double) n_unpaired / (n_unpaired + n_used + n_pcr_dupl)) << "); kept "
             << n_used + n_pcr_dupl << " read pairs in "
             << n_remaining_loci << " loci.\n";
        if (rm_pcr_duplicates) {
            cout << "Removed " << n_pcr_dupl
                 << " read pairs whose insert length had already been seen in the same sample as putative PCR duplicates ("
                 << as_percentage((double) n_pcr_dupl / (n_pcr_dupl + n_used))
                 << "); kept " << n_used << " read pairs.\n";
        } else {
            assert(n_pcr_dupl == 0);
        }
    } else {
        assert(gt_stats.n_unpaired_reads_rm() == 0 && gt_stats.n_read_pairs_pcr_dupl() == 0);
    }

    // Report statistics on genotyping and haplotyping.
    size_t n_hap_attempts = hap_stats.n_halplotyping_attempts();
    size_t n_hap_pairs = hap_stats.n_consistent_hap_pairs();
    OnlineMeanVar ecov_mean;
    double ecov_min = DBL_MAX;
    double ecov_max = 0.0;
    for (const GenotypeStats::PerSampleStats& sample : gt_stats.per_sample_stats) {
        double cov = (double) sample.ns_weighted_n_read_pairs_used / sample.ns_cumsum;
        ecov_mean.increment(cov);
        if (cov < ecov_min)
            ecov_min = cov;
        if (cov > ecov_max)
            ecov_max = cov;
    }
    o_fp1 << "\n"
         << "Genotyped " << gt_stats.n_genotyped_loci << " loci:\n";
    if (gt_stats.n_genotyped_loci == 0) {
        cerr << "Error: There wasn't any locus to genotype; check input/arguments.\n";
        throw exception();
    }
    o_fp1 << "  effective per-sample coverage: mean=" << ecov_mean.mean() << "x, stdev="
         << ecov_mean.sd_p() << "x, min=" << ecov_min << "x, max=" << ecov_max << "x\n"
         << "  mean number of sites per locus: " << gt_stats.mean_n_sites_per_loc() << "\n"
         << "  a consistent phasing was found for " << n_hap_pairs << " of out " << n_hap_attempts
         << " (" << as_percentage((double) n_hap_pairs / n_hap_attempts)
         << ") diploid loci needing phasing\n";

    // effective_coverages_per_sample
    logger->x << "\n"
                 "BEGIN effective_coverages_per_sample\n"
                 "# For mean_cov_ns, the coverage at each locus is weighted by the number of\n"
                 "# samples present at that locus (i.e. coverage at shared loci counts more).\n";
    logger->x << "sample"
                 "\tn_loci"
                 "\tn_used_fw_reads"
                 "\tmean_cov"
                 "\tmean_cov_ns";
    if (rm_unpaired_reads) {
        logger->x << "\tn_unpaired_reads";
        if (rm_pcr_duplicates) {
            logger->x << "\tn_pcr_dupl_pairs"
                         "\tpcr_dupl_rate";
        }
    }
    logger->x << '\n';
    assert(gt_stats.per_sample_stats.size() == bam_mpopi->samples().size());
    for (size_t sample=0; sample<bam_mpopi->samples().size(); ++sample) {
        const GenotypeStats::PerSampleStats& stats = gt_stats.per_sample_stats[sample];
        x_fp3 << bam_mpopi->samples()[sample].name
           << '\t' << stats.n_loci_with_sample
           << '\t' << stats.n_read_pairs_used
           << '\t' << (double) stats.n_read_pairs_used / stats.n_loci_with_sample
           << '\t' << (double) stats.ns_weighted_n_read_pairs_used / stats.ns_cumsum;
        if (rm_unpaired_reads) {
            x_fp3 << '\t' << stats.n_unpaired_reads;
            if (rm_pcr_duplicates) {
                x_fp3 << '\t' << stats.n_read_pairs_pcr_dupl
                      << '\t' << (double) stats.n_read_pairs_pcr_dupl / (stats.n_read_pairs_pcr_dupl + stats.n_read_pairs_used);
            }
        }
        x_fp3 << '\n';
    }
    logger->x << "END effective_coverages_per_sample\n";

    // pcr_clone_size_distrib
    if(rm_pcr_duplicates) {
        const vector<size_t>& cz_d = gt_stats.pcr_clone_size_distrib;
        assert(!cz_d.empty());
        assert(cz_d[0] == 0); // (We can't observe a clone of 0 reads.)
        logger->x << "\n"
                  << "BEGIN pcr_clone_size_distrib\n"
                  <<"clone_size\tn_clones\tn_reads\n";
        for(size_t clone_size=1; clone_size<cz_d.size(); ++clone_size)
            logger->x << clone_size
                      << '\t' << cz_d[clone_size]
                      << '\t' << cz_d[clone_size] * clone_size
                      << '\n';
        logger->x << "END pcr_clone_size_distrib\n";
    }

    // phasing_rates_samples
    logger->x << "\n"
              << "BEGIN phasing_rates_per_sample\n"
              << "sample\tn_gts\tn_multisnp_hets\tn_phased\tmisphasing_rate"
              #ifdef DEBUG
              << "\tn_phased_2ndpass"
              #endif
              << "\n";

    assert(hap_stats.per_sample_stats.size() == bam_mpopi->samples().size());
    for (size_t sample=0; sample<bam_mpopi->samples().size(); ++sample) {
        const HaplotypeStats::PerSampleStats& stats = hap_stats.per_sample_stats[sample];

        x_fp3 << bam_mpopi->samples()[sample].name
           << '\t' << stats.n_diploid_loci
           << '\t' << stats.n_hets_2snps
           << '\t' << stats.n_phased
           << '\t' << 1.0 - (double) stats.n_phased / stats.n_hets_2snps
           #ifdef DEBUG
           << '\t' << stats.n_phased_2ndpass
           #endif
           << '\n';
    }
    logger->x << "END phasing_rates_per_sample\n";

    // phasing_rates_loci
    #ifdef DEBUG
    logger->x << "\n"
              << "BEGIN phasing_rates_loci\n"
              << "n_hets_2snps\tn_bad_hets\tn_loci\n";
    for (auto& elem : hap_stats.n_badly_phased_samples)
        logger->x << elem.first.second
                  << '\t' << elem.first.second - elem.first.first
                  << '\t' << elem.second
                  << '\n';
    logger->x << "END phasing_rates_loci\n";
    #endif

    // Report clockings.
    {
        double ll  = t_parallel.elapsed();
        double v   = t_writing_vcf.elapsed();

        double r   = t_threads_totals.reading.elapsed() / num_threads;
        double p   = t_threads_totals.processing.elapsed() / num_threads;
        double w_f = t_threads_totals.writing_fa.elapsed() / num_threads;
        double w_v = t_threads_totals.writing_vcf.elapsed() / num_threads;
        double w_d = t_threads_totals.writing_details.elapsed() / num_threads;
        double w_b = t_threads_totals.writing_bams.elapsed() / num_threads;

        double ppr = t_threads_totals.processing_pre_alns.elapsed() / num_threads;
        double rn   = t_threads_totals.rm_Ns.elapsed() / num_threads;
        double as   = t_threads_totals.assembling.elapsed() / num_threads;
        double ia   = t_threads_totals.init_alignments.elapsed() / num_threads;
        double al   = t_threads_totals.aligning.elapsed() / num_threads;
        double me   = t_threads_totals.merge_paired_reads.elapsed() / num_threads;
        double ppo = t_threads_totals.processing_post_alns.elapsed() / num_threads;
        double rr = t_threads_totals.rm_reads.elapsed() / num_threads;
        double cnt = t_threads_totals.counting_nts.elapsed() / num_threads;
        double g   = t_threads_totals.genotyping.elapsed() / num_threads;
        double h   = t_threads_totals.haplotyping.elapsed() / num_threads;
        double u   = t_threads_totals.cpt_consensus.elapsed() / num_threads;
        double b_v = t_threads_totals.building_vcf.elapsed() / num_threads;
        double b_f = t_threads_totals.building_fa.elapsed() / num_threads;
        double b_b = t_threads_totals.building_bams.elapsed() / num_threads;

        double c = t_parallel.consumed()
                 + t_threads_totals.reading.consumed() / num_threads
                 + t_threads_totals.processing.consumed() / num_threads
                 + t_threads_totals.writing_fa.consumed() / num_threads
                 + t_threads_totals.writing_vcf.consumed() / num_threads
                 + t_threads_totals.writing_details.consumed() / num_threads
                 + t_threads_totals.writing_bams.consumed() / num_threads
                 + t_threads_totals.processing_pre_alns.consumed() / num_threads
                 + t_threads_totals.rm_Ns.consumed() / num_threads
                 + t_threads_totals.assembling.consumed() / num_threads
                 + t_threads_totals.init_alignments.consumed() / num_threads
                 + t_threads_totals.aligning.consumed() / num_threads
                 + t_threads_totals.merge_paired_reads.consumed() / num_threads
                 + t_threads_totals.processing_post_alns.consumed() / num_threads
                 + t_threads_totals.rm_reads.consumed() / num_threads
                 + t_threads_totals.counting_nts.consumed() / num_threads
                 + t_threads_totals.genotyping.consumed() / num_threads
                 + t_threads_totals.haplotyping.consumed() / num_threads
                 + t_threads_totals.cpt_consensus.consumed() / num_threads
                 + t_threads_totals.building_vcf.consumed() / num_threads
                 + t_threads_totals.building_fa.consumed() / num_threads
                 + t_threads_totals.building_bams.consumed() / num_threads
                 + t_writing_vcf.consumed()
                 ;

        x_fp1 << "\n"
           << "BEGIN clockings\n"
           << "Num. threads: " << num_threads << "\n"
           << "Parallel time: " << ll << "\n"
           << "Average thread time spent:\n"
           << std::setw(8) << r  << "  reading (" << as_percentage(r / ll) << ")\n"
           << std::setw(8) << p << "  processing (" << as_percentage(p / ll) << ")\n";
        if (as != 0.0)
            // De novo mode & paired-ends.
            x_fp1
               << std::setw(16) << ppr << " pre-alignments block (" << as_percentage(ppr / ll) << ")\n"
               << std::setw(16) << rn << "  reformatting fw-reads (" << as_percentage(rn / ll) << ")\n"
               << std::setw(16) << as << "  assembling (" << as_percentage(as / ll) << ")\n"
               << std::setw(16) << ia << "  initializing alignments (" << as_percentage(ia / ll) << ")\n"
               << std::setw(16) << al << "  aligning (" << as_percentage(al / ll) << ")\n"
               << std::setw(16) << me << "  merging read pairs (" << as_percentage(me / ll) << ")\n"
               ;
        x_fp1
           << std::setw(16) << ppo << " post-alignments block (" << as_percentage(ppo / ll) << ")\n"
           << std::setw(16) << rr << "  filtering reads (" << as_percentage(rr / ll) << ")\n"
           << std::setw(16) << cnt << "  counting nucleotides (" << as_percentage(cnt / ll) << ")\n"
           << std::setw(16) << g << "  genotyping (" << as_percentage(g / ll) << ")\n"
           << std::setw(16) << h << "  haplotyping (" << as_percentage(h / ll) << ")\n"
           << std::setw(16) << u << "  computing consensus (" << as_percentage(u / ll) << ")\n"
           << std::setw(16) << b_f << "  building_fa (" << as_percentage(b_f / ll) << ")\n"
           << std::setw(16) << b_v << "  building_vcf (" << as_percentage(b_v / ll) << ")\n";
        if (bam_output)
            x_fp1 << std::setw(16) << b_b << "  building_bam (" << as_percentage(b_b / ll) << ")\n"
                << std::setw(8) << w_b << "  writing_bam (" << as_percentage(w_b / ll) << ")\n";
        x_fp1 << std::setw(8) << w_f << "  writing_fa (" << as_percentage(w_f / ll) << ")\n"
           << std::setw(8) << w_v << "  writing_vcf (" << as_percentage(w_v / ll) << ")\n";
        if (detailed_output)
            x_fp1 << std::setw(8) << w_d << "  writing_details (" << as_percentage(w_d / ll) << ")\n";
        x_fp1 << std::setw(8) << c << "  clocking (" << as_percentage(c / ll) << ")\n"
           << "Total time spent writing vcf: " << v << " (" << as_percentage(v / ll) << ")\n"
           << "VCFwrite block size: mean=" << (double) gt_stats.n_genotyped_loci / n_writes
               << "(n=" << n_writes << "); max=" << max_size_before_write << "\n"
           << "END clockings\n";
    }

    #ifdef DEBUG
    const MarukiLowModel* m = dynamic_cast<const MarukiLowModel*>(model.get());
    if (m != NULL) {
        // Report how often the "underflow" likelihood equations were used
        cout << "\n"
             << "DEBUG: marukilow: calc_ln_weighted_sums: "
             << m->n_wsum_tot() << " calls, "
             << m->n_wsum_underflows() << " underflows ("
             << as_percentage(m->n_wsum_underflows(), m->n_wsum_tot())
             << ").\n";
        cout << "DEBUG: marukilow: mean_err_rate: "
             << m->mean_err_rate() << '\n';
    }
    #endif

    //
    // Cleanup & return.
    //

    if(bgzf_close(o_gzfasta_f) != 0)
        throw ios::failure("bgzf_close");
    o_vcf_f->file().close();
    o_vcf_f.reset();
    for (unique_ptr<Bam>& b : bam_of_ptrs) {
        b->close();
        b.reset();
    }
    if (o_details_f) {
        o_details_f->close();
        o_details_f.reset();
    }
    model.reset();
    if (dbg_write_hapgraphs)
        o_hapgraphs_f << "}\n";
    cout << "\ngstacks is done.\n";
    return 0;

} catch (exception& e) {
    return stacks_handle_exceptions(e);
}
}


//
// SnpAlleleCooccurrenceCounter
// ============================
//

const size_t& SnpAlleleCooccurrenceCounter::at(size_t snp_i1, Nt2 snp1_allele, size_t snp_i2, Nt2 snp2_allele) const {
    assert(snp_i1 < snp_i2);
    return cooccurrences_[snp_i1*n_snps_+snp_i2][size_t(snp1_allele)][size_t(snp2_allele)];
}

void SnpAlleleCooccurrenceCounter::clear() {
    for(size_t i=0; i<n_snps_; ++i)
        for(size_t j=i+1; j<n_snps_; ++j)
            cooccurrences_[i*n_snps_+j] = array<array<size_t,4>,4>();
}

//
// GenotypeStats, HaplotypeStats & ContigStats
// ===============
//

GenotypeStats& GenotypeStats::operator+= (const GenotypeStats& other) {
    this->n_genotyped_loci += other.n_genotyped_loci;
    this->n_sites_tot      += other.n_sites_tot;
    if (pcr_clone_size_distrib.size() < other.pcr_clone_size_distrib.size())
        pcr_clone_size_distrib.resize(other.pcr_clone_size_distrib.size());
    for (size_t i=0; i<other.pcr_clone_size_distrib.size(); ++i)
        pcr_clone_size_distrib[i] += other.pcr_clone_size_distrib[i];
    for (size_t sample=0; sample<per_sample_stats.size(); ++sample) {
        per_sample_stats[sample].n_unpaired_reads += other.per_sample_stats[sample].n_unpaired_reads;
        per_sample_stats[sample].n_read_pairs_pcr_dupl += other.per_sample_stats[sample].n_read_pairs_pcr_dupl;
        per_sample_stats[sample].n_read_pairs_used += other.per_sample_stats[sample].n_read_pairs_used;
        per_sample_stats[sample].n_loci_with_sample += other.per_sample_stats[sample].n_loci_with_sample;
        per_sample_stats[sample].ns_cumsum += other.per_sample_stats[sample].ns_cumsum;
        per_sample_stats[sample].ns_weighted_n_read_pairs_used += other.per_sample_stats[sample].ns_weighted_n_read_pairs_used;
    }
    return *this;
}

HaplotypeStats& HaplotypeStats::operator+= (const HaplotypeStats& other) {
    for (auto elem : other.n_badly_phased_samples)
        this->n_badly_phased_samples[elem.first] += elem.second;
    assert(per_sample_stats.size() == other.per_sample_stats.size());
    for (size_t sample=0; sample<per_sample_stats.size(); ++sample) {
        per_sample_stats[sample].n_diploid_loci += other.per_sample_stats[sample].n_diploid_loci;
        per_sample_stats[sample].n_hets_2snps += other.per_sample_stats[sample].n_hets_2snps;
        per_sample_stats[sample].n_phased += other.per_sample_stats[sample].n_phased;
        per_sample_stats[sample].n_phased_2ndpass += other.per_sample_stats[sample].n_phased_2ndpass;
    }
    return *this;
}

size_t HaplotypeStats::n_halplotyping_attempts() const {
    size_t n1 = 0;
    for (auto& s: per_sample_stats)
        n1 += s.n_hets_2snps;
    size_t n2 = 0;
    for (auto& e: n_badly_phased_samples)
        n2+= e.second * e.first.second;
    assert(n1==n2);
    return n1;
}

size_t HaplotypeStats::n_consistent_hap_pairs() const {
    size_t n1 = 0;
    for (auto& s: per_sample_stats)
        n1 += s.n_phased;
    size_t n2 = 0;
    for (auto& e: n_badly_phased_samples)
        n2 += e.second * e.first.first;
    assert(n1==n2);
    return n1;
}

ContigStats& ContigStats::operator+= (const ContigStats& other) {
    this->n_nonempty_loci           += other.n_nonempty_loci;
    this->n_loci_w_pe_reads         += other.n_loci_w_pe_reads;
    this->n_loci_almost_no_pe_reads += other.n_loci_almost_no_pe_reads;
    this->n_loci_pe_graph_not_dag   += other.n_loci_pe_graph_not_dag;
    this->n_loci_pe_graph_fixed_to_dag += other.n_loci_pe_graph_fixed_to_dag;
    this->length_ctg_tot            += other.length_ctg_tot;
    this->n_aln_reads               += other.n_aln_reads;
    this->n_tot_reads               += other.n_tot_reads;
    this->n_overlaps                += other.n_overlaps;
    this->length_overlap_tot        += other.length_overlap_tot;
    this->length_olapd_loci_tot     += other.length_olapd_loci_tot;
    this->insert_length_olap_mv    += other.insert_length_olap_mv;

    return *this;
}

//
// LocData
// =======
//
void LocData::clear() {
    id = -1;
    pos.clear();
    mpopi = NULL;
    ctg_status = unknown;
    olap_length = SIZE_MAX;
    o_vcf.clear();
    o_fa.clear();
    o_details.clear();
    details_ss.clear();
    details_ss.str(string());
}

//
// LocusProcessor
// ==============
//

void
LocusProcessor::process(CLocReadSet& loc)
{
    timers_.processing_pre_alns.restart();
    loc_.clear();
    if (dbg_denovo_min_loc_samples > 0
            && loc.n_samples() < dbg_denovo_min_loc_samples)
        loc.clear();
    if (loc.reads().empty()) {
        timers_.processing_pre_alns.update();
        return;
    }
    if (detailed_output)
        loc_.details_ss << "BEGIN locus " << loc.id() << "\n";

    ++ctg_stats_.n_nonempty_loci;
    if (!loc.pe_reads().empty())
        ++ctg_stats_.n_loci_w_pe_reads;
    loc_.id = loc.id();
    loc_.pos = loc.pos();
    loc_.mpopi = &loc.mpopi();

    CLocAlnSet aln_loc;
    aln_loc.reinit(loc_.id, loc_.pos, loc_.mpopi);

    //
    // Remove N's in the forward reads.
    //
    timers_.rm_Ns.restart();
    for (SRead& r : loc.reads())
        r.seq.remove_Ns();
    timers_.rm_Ns.update();

    //
    // Assemble a contig.
    //
    timers_.assembling.restart();
    if (detailed_output)
        loc_.details_ss << "BEGIN contig\n";
    DNASeq4 ctg = assemble_locus_contig(loc.reads(), loc.pe_reads());
    timers_.assembling.update();
    timers_.processing_pre_alns.update();
    if (ctg.empty()) {
        if (detailed_output) {
            loc_.details_ss << "contig_assembly_failed\n"
                << "END contig\n"
                << "END locus " << loc_.id << '\n';
            loc_.o_details = loc_.details_ss.str();
        }
        return;
    }
    aln_loc.ref(move(ctg));
    timers_.assembling.update();

    timers_.init_alignments.restart();
    SuffixTree* stree = new SuffixTree(aln_loc.ref());
    stree->build_tree();
    GappedAln aligner;
    AlignRes aln_res;
    if (!loc.pe_reads().empty()) {
        //
        // Check that the contig starts at the cutsite (and not in the
        // adapter), by aligning one FW read.
        //
        for (SRead& r : loc.reads()) {
            if (!align_reads_to_contig(stree, &aligner, r.seq, aln_res))
                continue;
            // Read did align. Check start position & break.
            if (aln_res.subj_pos > 0) {
                DNASeq4 new_ctg;
                DNASeq4::iterator itr = aln_loc.ref().begin();
                for (size_t i=0; i<aln_res.subj_pos; ++i) {
                    assert(itr != aln_loc.ref().end());
                    ++itr;
                }
                new_ctg.append(itr, aln_loc.ref().end());
                if (detailed_output)
                    loc_.details_ss
                        << "orig_contig\t" << aln_loc.ref() << '\n'
                        << "first_fw_read_aln\t" << r.name << "\tpos:" << aln_res.subj_pos << '\n';
                aln_loc.ref(move(new_ctg));
                delete stree;
                stree = new SuffixTree(aln_loc.ref());
                stree->build_tree();
            }
            break;
        }
    }
    ctg_stats_.length_ctg_tot += aln_loc.ref().length();
    if (detailed_output)
        loc_.details_ss << "contig\t" << aln_loc.ref() << '\n'
            << "END contig\n";
    timers_.init_alignments.update();

    //
    // Align the reads.
    //
    timers_.aligning.restart();
    this->ctg_stats_.n_tot_reads += loc.reads().size();
    this->ctg_stats_.n_tot_reads += loc.pe_reads().size();
    if (detailed_output)
        loc_.details_ss << "BEGIN pe_alns\n";
    for (SRead& r : loc.reads()) { // FORWARD READS
        if(add_read_to_aln(aln_loc, aln_res, move(r), &aligner, stree)) {
            this->ctg_stats_.n_aln_reads++;
            if (detailed_output)
                loc_.details_ss << "fw_aln_local"
                                << '\t' << aln_loc.reads().back().name
                                << '\t' << aln_res.subj_pos + 1 << ':' << aln_res.cigar
                                << '\n';
        } else {
            if (detailed_output)
                // `r` hasn't been moved.
                loc_.details_ss << "fw_aln_fail\t" << r.name << '\n';
        }
    }
    for (SRead& r : loc.pe_reads()) {
        if (add_read_to_aln(aln_loc, aln_res, move(r), &aligner, stree)) {
            this->ctg_stats_.n_aln_reads++;
            if (loc_.ctg_status == LocData::overlapped) {
                // Record the insert length. (Insert lengths are just based on
                // reverse reads but they should all have matching foward reads
                // starting at the restriction site so it makes sense.)
                const Cigar& cigar = aln_loc.reads().back().cigar;
                size_t insert_length = aln_loc.ref().length();
                assert(!cigar.empty());
                if (cigar.back().first == 'D') {
                    insert_length -= cigar.back().second;
                    if (cigar.size() >=2 && (*++cigar.rbegin()).first == 'I')
                        // Soft clipping on the far end.
                        insert_length += (*++cigar.rbegin()).second;
                }
                ctg_stats_.insert_length_olap_mv.increment(insert_length);
            }
            if (detailed_output)
                loc_.details_ss << "pe_aln_local"
                                << '\t' << aln_loc.reads().back().name
                                << '\t' << aln_res.subj_pos + 1 << ':' << aln_res.cigar
                                << '\n';
        } else {
            if (detailed_output)
                // `r` hasn't been moved.
                loc_.details_ss << "pe_aln_fail\t" << r.name << '\n';
        }
    }
    if (detailed_output)
        loc_.details_ss << "END pe_alns\n";
    delete stree;
    timers_.aligning.update();

    if (bam_output) {
        timers_.building_bams.restart();
        assert(input_type == GStacksInputT::denovo_popmap || input_type == GStacksInputT::denovo_merger);
        const bool merged_input = input_type == GStacksInputT::denovo_merger;
        if (loc_.o_bam.empty())
            loc_.o_bam.resize(merged_input ? 1 : loc_.mpopi->n_samples());
        aln_loc.sort_by_alignment_offset();
        for (const SAlnRead& read : aln_loc.reads()) {
            BamRecord rec;
            assert(read.name.length() >= 2 && *++read.name.rbegin() == '/');
            Cigar c = read.cigar;
            size_t offset = 0;
            assert(!c.empty());
            if (c.back().first == 'D')
                c.pop_back();
            else if (c.back().first == 'I' && c.size() >= 2 && (*----c.end()).first == 'D')
                c.erase(----c.end());
            if (c[0].first == 'D') {
                offset = c[0].second;
                c.erase(c.begin());
            } else if (c[0].first == 'I' && c.size() >= 2 && c[1].first == 'D') {
                offset = c[1].second;
                c.erase(++c.begin());
            }
            assert(bam_cloc_reader_);
            rec.assign(
                read.name.substr(0, read.name.length() - 2),
                *read.name.rbegin() == '1' ? BAM_FREAD1 : BAM_FREAD2,
                bam_cloc_reader_->id2target(aln_loc.id()),
                offset,
                c,
                read.seq,
                read.sample
                );
            loc_.o_bam.at(merged_input ? 0 : read.sample).push_back(move(rec));
        }
        timers_.building_bams.update();
    }

    timers_.merge_paired_reads.restart();
    aln_loc.merge_paired_reads();
    timers_.merge_paired_reads.update();
    loc.clear();
    timers_.processing_pre_alns.update();
    if (aln_loc.reads().empty()) {
       if (detailed_output)
           loc_.details_ss << "no_reads_were_aligned\n";
        return;
    }
    process(aln_loc);
}

void
LocusProcessor::process(CLocAlnSet& aln_loc)
{
    assert(!aln_loc.reads().empty());
    timers_.processing_post_alns.restart();
    if (input_type == GStacksInputT::denovo_popmap || input_type == GStacksInputT::denovo_merger) {
        // Called from process(CLocReadSet&).
        assert(this->loc_.id == aln_loc.id());
    } else {
        this->loc_.clear();
        this->loc_.id = aln_loc.id();
        this->loc_.pos = aln_loc.pos();
        this->loc_.mpopi = &aln_loc.mpopi();
        if (detailed_output)
            this->loc_.details_ss << "BEGIN locus " << loc_.id << "\n";
    }

    //
    // Remove unpaired reads/PCR duplicates.
    //
    timers_.rm_reads.restart();
    if (rm_unpaired_reads) {
        for (size_t sample=0; sample<gt_stats_.per_sample_stats.size(); ++sample)
            gt_stats_.per_sample_stats[sample].n_unpaired_reads += aln_loc.sample_reads(sample).size();
        aln_loc.remove_unmerged_reads(detailed_output ? &this->loc_.details_ss : NULL);
        for (size_t sample=0; sample<gt_stats_.per_sample_stats.size(); ++sample)
            gt_stats_.per_sample_stats[sample].n_unpaired_reads -= aln_loc.sample_reads(sample).size();
        if(aln_loc.reads().empty()) {
            loc_.clear();
            timers_.processing_post_alns.update();
            return;
        }
    }
    if (rm_pcr_duplicates) {
        assert(rm_unpaired_reads);
        for (size_t sample=0; sample<gt_stats_.per_sample_stats.size(); ++sample)
            gt_stats_.per_sample_stats[sample].n_read_pairs_pcr_dupl += aln_loc.sample_reads(sample).size();
        aln_loc.remove_pcr_duplicates(&gt_stats_.pcr_clone_size_distrib,
                                      detailed_output ? &this->loc_.details_ss : NULL);
        for (size_t sample=0; sample<gt_stats_.per_sample_stats.size(); ++sample)
            gt_stats_.per_sample_stats[sample].n_read_pairs_pcr_dupl -= aln_loc.sample_reads(sample).size();
    }
    timers_.rm_reads.update();
    assert(!aln_loc.reads().empty());
    ++gt_stats_.n_genotyped_loci;
    size_t n_samples = aln_loc.n_samples();
    for (size_t sample=0; sample<gt_stats_.per_sample_stats.size(); ++sample) {
        if (!aln_loc.sample_reads(sample).empty()) {
            auto& s = gt_stats_.per_sample_stats[sample];
            ++s.n_loci_with_sample;
            s.n_read_pairs_used += aln_loc.sample_reads(sample).size();
            s.ns_cumsum += n_samples;
            s.ns_weighted_n_read_pairs_used += n_samples * aln_loc.sample_reads(sample).size();
        }
    }
    if (detailed_output) {
        loc_.details_ss << "BEGIN aln_matrix\n";
        for (const SAlnRead& read : aln_loc.reads())
            loc_.details_ss << read.name << '\t' << loc_.mpopi->samples()[read.sample].name
                            << '\t' << read.cigar << '\n';
        loc_.details_ss << "END aln_matrix\n";
    }

    //
    // Get the nucleotide counts at each position.
    //
    timers_.counting_nts.restart();
    vector<SiteCounts> depths;
    depths.reserve(aln_loc.ref().length());
    for (CLocAlnSet::site_iterator site (aln_loc); bool(site); ++site) {
        depths.push_back(site.counts());
        if (depths.back().tot.sum() > 0)
            ++gt_stats_.n_sites_tot; // Sites "with data".
    }
    assert(depths.size() == aln_loc.ref().length());
    timers_.counting_nts.update();

    //
    // Call SNPs.
    //
    timers_.genotyping.restart();
    vector<SiteCall> calls;
    calls.reserve(aln_loc.ref().length());
    for (const SiteCounts& site_depths : depths) {
        calls.push_back(model->call(site_depths));
        // Make sure our SNP and genotype calls are consistent (discard SNPs with
        // a MAC of 0, or without calls at all).
        if (calls.back().alleles().size() > 1)
            calls.back().filter_mac(1);
    }
    timers_.genotyping.update();

    //
    // Call haplotypes; amend SNP/genotype calls.
    //
    vector<map<size_t,PhasedHet>> phase_data;
    if (!dbg_no_haplotypes) {
        timers_.haplotyping.restart();
        phase_hets(phase_data, calls, aln_loc, hap_stats_);
        timers_.haplotyping.update();
    }

    //
    // Determine the consensus sequence.
    //
    timers_.cpt_consensus.restart();
    DNASeq4 consensus (aln_loc.ref().length());
    assert(calls.size() == aln_loc.ref().length());
    assert(depths.size() == aln_loc.ref().length());
    for (size_t i = 0; i < aln_loc.ref().length(); ++i) {
        if (!calls[i].alleles().empty()) {
            consensus.set(i, calls[i].most_frequent_allele());
        } else {
            // Use the majority-rule nucleotide.
            // (For the high/low Maruki" models this actually only happens when
            // there is no coverage; for the Hohenlohe model it may also happen
            // when there isn't a single significant call.)
            pair<size_t,Nt2> nt = depths[i].tot.sorted()[0];
            if (nt.first > 0)
                consensus.set(i, nt.second);
        }
    }
    aln_loc.ref(move(consensus));
    timers_.cpt_consensus.update();

    //
    // Create the outputs.
    //
    write_one_locus(aln_loc, depths, calls, phase_data);

    if (detailed_output) {
        loc_.details_ss << "END locus " << loc_.id << "\n";
        loc_.o_details = loc_.details_ss.str();
    }
    if (dbg_write_alns) {
        #pragma omp critical(dbg_alns)
        o_aln_f << "BEGIN " << loc_.id << "\n"
                << aln_loc
                << "\nEND " << loc_.id << "\n";
    }
    timers_.processing_post_alns.update();
}

bool
LocusProcessor::align_reads_to_contig(SuffixTree *st, GappedAln *g_aln, DNASeq4 &enc_query, AlignRes &aln_res) const
{
    vector<STAln> alns, final_alns;
    vector<pair<size_t, size_t> > step_alns;
    string      query  = enc_query.str();
    const char *q      = query.c_str();
    const char *q_stop = q + query.length();
    size_t      q_pos  = 0;
    size_t      id     = 0;
    char c[id_len];

    do {
        step_alns.clear();

        q_pos = q - query.c_str();

        st->align(q, step_alns);

        if (step_alns.size() == 0 || step_alns.size() > max_fragment_alns) {
            q++;
        } else {
            for (uint i = 0; i < step_alns.size(); i++)
                alns.push_back(STAln(id, q_pos, step_alns[i].first, step_alns[i].second));
            q += step_alns[0].second + 1;
            id++;
        }
    } while (q < q_stop);

    //
    // No alignments to the suffix tree were found.
    //
    if (alns.size() == 0)
        return false;

    //
    // Perfect alignmnet to the suffix tree. Return result.
    //
    if (alns.size() == 1 && alns[0].aln_len == query.length()) {
        snprintf(c, id_len, "%luM", query.length());
        aln_res.cigar    = c;
        aln_res.subj_pos = alns[0].subj_pos;
        return true;
    }

    //
    // Find a consistent set of suffix tree alignments, ordered in a directed, acyclic grapgh (DAG).
    //
    this->suffix_tree_hits_to_dag(query.length(), alns, final_alns);

    g_aln->init(query.length(), st->seq_len(), true);
    if (g_aln->align_constrained(query, st->seq_str(), final_alns)) {
        aln_res = g_aln->result();
    }

    if (aln_res.pct_id < min_aln_cov)
        return false;

    return true;
}

int
LocusProcessor::suffix_tree_hits_to_dag(size_t query_len, vector<STAln> &alns, vector<STAln> &final_alns) const
{
    //
    // 1. Sort the alignment fragments so they are primarily ordered by subject position, secondarily by query position.
    //
    sort(alns.begin(), alns.end(), [](STAln a, STAln b)
         {
            if (a.subj_pos == b.subj_pos)
                return a.query_pos < b.query_pos;
            else
                return a.subj_pos < b.subj_pos;
         });

    int    gap_len, end_pos;
    double score, scale;
    bool   term;
    vector<size_t> term_nodes;
    //
    // 2. Traverse the list of fragments and add links to reachable nodes, recording the score.
    //    2.1. A successor node is reachable from the predecessor node if both the query and subject
    //         positions are advanced relative to the predecessor, and there is no sequence overlap
    //         between the two nodes.
    //    2.2. We weight the score by the number of 'gap' nucleotides between the two fragments to
    //         penalize fragments that are further away.
    //    2.3  Record the maximal link in the graph.
    //    2.4. Determine if a node is a terminal node.
    //
    for (uint i = 0; i < alns.size(); i++) {
        end_pos = alns[i].subj_pos + alns[i].aln_len - 1;
        term    = true;

        assert(end_pos > 0);

        for (uint j = i + 1; j < alns.size(); j++) {
            if (alns[j].query_pos > alns[i].query_pos &&
                alns[j].subj_pos  > alns[i].subj_pos  &&
                alns[j].subj_pos  > (uint) end_pos) {

                term    = false;
                gap_len = alns[j].subj_pos - end_pos - 1;
                scale   = gap_len > (int) query_len ? 1 : alns[i].aln_len * ((double) gap_len / (double) query_len);
                score   = alns[i].aln_len - scale; // Raw score.
                score  += alns[i].max._score;      // Score adjusted for the highest incoming node.

                assert(gap_len >= 0);
                assert(score   >= 0);

                alns[j].links.push_back(STLink(i, score));
                if (alns[j].max._index == j || score > alns[j].max._score)
                    alns[j].max = STLink(i, score);
            }
        }
        if (term)
            term_nodes.push_back(i);
    }

    //
    // 3. Find the terminal node with the highest score.
    //
    double max_score = alns[term_nodes.front()].max._score;
    size_t max_index = term_nodes.front();
    for (uint i = 1; i < term_nodes.size(); i++) {
        if (alns[term_nodes[i]].max._score > max_score) {
            max_index = term_nodes[i];
            max_score = alns[term_nodes[i]].max._score;
        }
    }

    //
    // 4. Backtrack to get the optimal path.
    //
    vector<size_t> optimal;

    if (max_score == 0 && max_index == 0) {
        //
        // None of the fragments were connected in the DAG, select the largest fragment.
        //
        max_score = alns[term_nodes.front()].aln_len;
        max_index = term_nodes.front();
        for (uint i = 1; i < term_nodes.size(); i++)
            if (alns[i].aln_len > max_score) {
                max_index = term_nodes[i];
                max_score = alns[term_nodes[i]].aln_len;
            }
        optimal.push_back(max_index);

    } else {

        uint n = max_index;
        while (alns[n].links.size() > 0) {
            optimal.push_back(n);
            n = alns[n].max._index;
        }
        optimal.push_back(n);
    }
    assert(optimal.size() > 0);

    final_alns.clear();
    for (int i = optimal.size() - 1; i >= 0; i--)
        final_alns.push_back(alns[optimal[i]]);

    return 0;
}

size_t
LocusProcessor::find_locus_overlap(SuffixTree *stree, GappedAln *g_aln, const DNASeq4 &se_consensus, string &overlap_cigar) const
{
    vector<STAln> alns, final_alns;
    vector<pair<size_t, size_t> > step_alns;
    AlignRes aln_res;

    string      query  = se_consensus.str();
    const char *q      = query.c_str();
    const char *q_stop = q + query.length();
    size_t      q_pos  = 0;
    size_t      id     = 1;

    do {
        step_alns.clear();

        q_pos = q - query.c_str();

        stree->align(q, step_alns);

        if (step_alns.size() == 0 || step_alns.size() > max_fragment_alns) {
            q++;
        } else {
            for (uint i = 0; i < step_alns.size(); i++)
                alns.push_back(STAln(id, q_pos, step_alns[i].first, step_alns[i].second));
            q += step_alns[0].second + 1;
            id++;
        }
    } while (q < q_stop);


    if (alns.size() == 0) {
        //
        // If no alignments have been found, search the tails of the query and subject for any overlap
        // that is too small to be picked up by the SuffixTree using the gapped aligner.
        //
        uint olap_bound = stree->min_aln() * 2;

        //
        // Create a gapped alignment, looking only at the tail of the two sequences, to determine the exact overlap.
        //
        g_aln->init(query.length(), stree->seq_len(), true);
        if (g_aln->align_region(query, stree->seq_str(),
                                query.length() - olap_bound, query.length() - 1, 0, olap_bound - 1)) {
            aln_res       = g_aln->result();
            overlap_cigar = aln_res.cigar;
        }
        // g_aln->dump_alignment(query, stree->seq_str());

    } else {
        size_t query_stop = alns.front().query_pos + alns.front().aln_len - 1;

        if ( alns.size() == 1 &&
            (query_stop == (query.length() - 1) && alns.front().subj_pos == 0)) {
            //
            // If a single alignment has been found, check for a perfect alignmnet to the suffix tree
            // that occupies the end of the query and the beginning of the subject.
            //
            char buf[id_len];
            snprintf(buf, id_len, "%luM", alns.front().aln_len);
            overlap_cigar  = buf;
            aln_res.pct_id = 1;

        } else {
            //
            // Otherwise, find a consistent set of suffix tree alignments, ordered in a directed, acyclic grapgh (DAG).
            //
            this->suffix_tree_hits_to_dag(query.length(), alns, final_alns);

            //
            // Create a gapped alignment, prefilling the suffix tree hits, to determine the exact overlap.
            //
            g_aln->init(query.length(), stree->seq_len(), true);
            if (g_aln->align_constrained(query, stree->seq_str(), final_alns)) {
                aln_res       = g_aln->result();
                overlap_cigar = aln_res.cigar;
            }
        }
    }

    Cigar cigar;
    parse_cigar(overlap_cigar.c_str(), cigar, true);

    if (cigar.size() == 0)
        return 0;

    // cerr << "Cigar: " << aln_res.cigar.c_str() << "; cigar_length_ref aka overlap: " << cigar_length_ref(cigar) << "\n";

    //
    // Assess the gapped alignment and make sure it is reasonable.
    //
    size_t overlap     = 0;
    size_t cigar_index = 0;
    double total_len   = 0.0;
    double align_len   = 0.0;

    if (cigar[cigar_index].first == 'S') {
        //
        // Here is an example where we want to ignore the initial softmasking repored for the SE contig,
        // but we must count the offset to the start of the PE contig as a set of mismatches.
        // CIGAR:  138S4M2S
        // SE CTG: TGCAGGAGATTTTCCTCTCTGCTGT...GGAGCTCGTTGGTGTTGGTGGCACTGCCCTGAGGTTTTCCTCATNNTCCCCAGGNCCCAGCGTCAGC
        // PE CTG:                                                                               GGCACAGATGTGTCACAGCTGC...CTTCTT
        //                                                                                                  ****
        cigar_index++;
        total_len += aln_res.subj_pos;

    } else {
        //
        // If the PE contig overlaps the entire SE contig, and there is sequence remaining on the PE contig, adjust the
        // overlap to include it, but it does not count as unaligned total length.
        // e.g.:
        // CIGAR:  68M1D6M
        // SE CTG:     TGCAGGCGCACCTGTAGGCCCCCGGACAGGAGGGGGTGGCTTCAAGCAGGGGCCAGCCCGAGGCCCCC-ACCCAC
        // PE CTG: TAGTTGCAGGCGCACCTGTAGGCCCCCGGACAGGAGGGGGTGGCTTCAAGCAGGGGCCAGCCCGAGGCCCCCCACCCACACCCAGCACACACTGGCC...
        //
        overlap += aln_res.subj_pos;
    }

    for (; cigar_index < cigar.size(); cigar_index++) {
        total_len += cigar[cigar_index].second;
        if (cigar[cigar_index].first == 'M')
            align_len += cigar[cigar_index].second;
    }

    //
    // Test that 80% of the overlapping sequence is aligned (includes gaps) AND
    // make sure that of those sequences aligned, 80% are identities.
    //
    if (align_len / total_len < overlap_min_pct_id ||
        aln_res.pct_id        < overlap_min_pct_id)
        return 0;

    overlap += align_len;

    //
    // Return overlap.
    //
    return overlap;
}

DNASeq4
LocusProcessor::assemble_locus_contig(
    const vector<SRead>& fw_reads,
    const vector<SRead>& pe_reads
) {
    //
    // List the sequences that will enter the graph.
    //
    array<vector<const DNASeq4*>, 3> seqs_to_assemble;
    array<const vector<SRead>*, 2> reads = {{&fw_reads, &pe_reads}};
    size_t fw = 0, pe = 1, both = 2;
    for (size_t i : {fw, pe}) {
        if (reads[i]->size() <= max_debruijn_reads) {
            for (const Read& r : *reads[i])
                seqs_to_assemble[i].push_back(&r.seq);
        } else {
            double seq_i = 0.0;
            double step = (double) reads[i]->size() / max_debruijn_reads;
            for (size_t j=0; j<max_debruijn_reads; ++j) {
                seqs_to_assemble[i].push_back(&(*reads[i])[(size_t)seq_i].seq);
                seq_i += step;
            }
        }
    }
    seqs_to_assemble[both] = seqs_to_assemble[fw];
    seqs_to_assemble[both].insert(
        seqs_to_assemble[both].end(),
        seqs_to_assemble[pe].begin(), seqs_to_assemble[pe].end());

    //
    // Build the graph.
    //
    Graph graph (km_length);
    graph.rebuild(seqs_to_assemble[both], min_km_count);
    if (graph.empty())
        return DNASeq4();
    if (dbg_write_gfa) {
        graph.dump_gfa(out_dir + "gstacks." + to_string(loc_.id) + ".spaths.gfa");
        graph.dump_gfa(out_dir + "gstacks." + to_string(loc_.id) + ".nodes.gfa", true);
    }

    //
    // Find the best path in the graph.
    //
    vector<const SPath*> best_path;
    do {
        if (graph.find_best_path(best_path))
            break;

        // The graph wasn't a DAG, try to fix it.
        if(graph.remove_cycles() && graph.find_best_path(best_path)) {
            ++ctg_stats_.n_loci_pe_graph_fixed_to_dag;
            if (detailed_output)
                loc_.details_ss << "fixed_to_dag\n";
            if (dbg_write_gfa) {
                graph.dump_gfa(out_dir + "gstacks." + to_string(loc_.id) + ".fixed.spaths.gfa");
                graph.dump_gfa(out_dir + "gstacks." + to_string(loc_.id) + ".fixed.nodes.gfa", true);
            }
            break;
        }

        // Contig contruction failed.
        ++ctg_stats_.n_loci_pe_graph_not_dag;
        if (detailed_output)
            loc_.details_ss << "not_dag\n";
        if (dbg_write_gfa_notdag) {
            graph.dump_gfa(out_dir + "gstacks." + to_string(loc_.id) + ".spaths.gfa");
            graph.dump_gfa(out_dir + "gstacks." + to_string(loc_.id) + ".nodes.gfa", true);
        }
        return DNASeq4();
    } while (false);
    if (detailed_output)
        loc_.details_ss << "is_dag\n";
    DNASeq4 contig = DNASeq4(SPath::contig_str(best_path.begin(), best_path.end(), km_length));

    //
    // Handle the case where there are paired-end reads and the main
    // corresponding subgraph is distinct from that of the forward reads.
    //
    if (!pe_reads.empty()) {
        //
        // Extract the components and map the fw & pe reads to them.
        //
        vector<vector<const SPath*>> components = graph.components();
        unordered_map<Kmer,size_t> component_kmers;
        for (size_t i=0; i<components.size(); ++i) {
            size_t n_nodes = 0;
            for (const SPath* sp : components[i]) {
                for (const Node* n=sp->first(); ; n=n->first_succ()) {
                    component_kmers[n->km()] = i;
                    ++n_nodes;
                    if (n == sp->last())
                        break;
                }
            }
        }
        array<vector<size_t>, 2> comp_depths;
        for (size_t i : {fw, pe}) {
            comp_depths[i] = vector<size_t>(components.size(), 0);
            for (const DNASeq4* s : seqs_to_assemble[i]) {
                Kmerizer kmers {km_length, *s};
                Kmer km;
                decltype(component_kmers.begin()) itr;
                while ((km = kmers.next()))
                    if ((itr = component_kmers.find(km)) != component_kmers.end())
                        ++comp_depths[i][itr->second];
            }
        }
        //
        // Check that there actually exist kmers for the forward and paired-end reads.
        //
        if (accumulate(comp_depths[fw].begin(), comp_depths[fw].end(), 0) == 0) {
            if (detailed_output)
                loc_.details_ss << "fw_reads_absent_from_graph\n";
            return DNASeq4();
        }
        if (accumulate(comp_depths[pe].begin(), comp_depths[pe].end(), 0) == 0) {
            if (detailed_output)
                loc_.details_ss << "pe_reads_absent_from_graph\n";
        } else {
            //
            // Check that the main forward and paired-end components are the same one.
            //
            array<size_t, 2> best_comp;
            for (size_t i : {fw, pe}) {
                size_t best_depth = 0;
                for (size_t c=0; c<components.size(); ++c) {
                    if (comp_depths[i][c] > best_depth) {
                        best_comp[i] = c;
                        best_depth = comp_depths[i][c];
                    }
                }
            }
            if (best_comp[fw] == best_comp[pe]) {
                loc_.ctg_status = LocData::overlapped;
                loc_.olap_length = km_length;
                ctg_stats_.n_overlaps++;
                ctg_stats_.length_overlap_tot += loc_.olap_length;
                ctg_stats_.length_olapd_loci_tot += contig.length();
                if (detailed_output)
                    loc_.details_ss << "one_debruijn_subgraph\n"
                        << "overlap\t" << km_length << '\n';
            } else {
                if (detailed_output)
                    loc_.details_ss << "two_debruijn_subgraphs\n";
                //
                // Determine the component the best path corresponds to.
                //
                Kmer km = Kmerizer(km_length, contig).next();
                assert(km);
                auto comp_itr = component_kmers.find(km);
                assert(comp_itr != component_kmers.end());
                size_t comp = comp_itr->second;
                //
                // Get two contigs, one for forward and one for paired-end.
                //
                DNASeq4 fw_contig, pe_contig;
                DNASeq4* todo_contig;
                const vector<const DNASeq4*>* todo_seqs;
                if (comp == best_comp[fw]) {
                    if (detailed_output)
                        loc_.details_ss << "best_path_is_best_fw_compo\n";
                    fw_contig = move(contig);
                    todo_contig = &pe_contig;
                    todo_seqs = &seqs_to_assemble[pe];
                } else if (comp == best_comp[pe]) {
                    if (detailed_output)
                        loc_.details_ss << "best_path_is_best_pe_compo\n";
                    pe_contig = move(contig);
                    todo_contig = &fw_contig;
                    todo_seqs = &seqs_to_assemble[fw];
                } else {
                    if (detailed_output)
                        loc_.details_ss << "best_path_is_neither_best_fw_compo_nor_best_pe_compo\n";
                    ++ctg_stats_.n_loci_pe_graph_not_dag;
                    return DNASeq4();
                }
                graph.rebuild(*todo_seqs, min_km_count);
                if (graph.empty()) {
                    // DOES_NOT_HAPPEN; // May fail, c.f. mailing list 2019-01-18.
                    // We can arrive here in some limit cases. Especially, if the forward and reverse
                    // regions overlap but the graph is still disconnected, which may happen e.g. if
                    // there are very few reads and some sequencing errors/SNPs (e.g. two 150bp forward
                    // reads and one 150bp reverse read with a 50bp overlap, and the two forward reads
                    // differ at position 100 and min_km_count is 2; with enough errors in the forward
                    // reads it's probably also possible to make this happen with the graph's best path
                    // being in the reverse region.)
                    loc_.details_ss << "limit_case_the_other_compo_disappeared\n";
                    return DNASeq4();
                }
                if (!graph.find_best_path(best_path)
                        && !(graph.remove_cycles() && graph.find_best_path(best_path)))
                    DOES_NOT_HAPPEN;
                *todo_contig = DNASeq4(SPath::contig_str(best_path.begin(), best_path.end(), km_length));
                if(detailed_output)
                    loc_.details_ss << "contig_fw\t" << fw_contig << '\n'
                                    << "contig_pe\t" << pe_contig << '\n';
                //
                // Try to overlap the two contigs.
                //
                SuffixTree stree {pe_contig};
                stree.build_tree();
                GappedAln aligner;
                AlignRes  aln_res;
                string overlap_cigar;
                int overlap;
                if (dbg_no_overlaps)
                    overlap = 0;
                else
                    overlap = this->find_locus_overlap(&stree, &aligner, fw_contig, overlap_cigar);
                if (overlap > 0) {
                    if(detailed_output)
                        loc_.details_ss << "overlap\t" << overlap << '\t' << overlap_cigar << "\n";
                    this->loc_.ctg_status = LocData::overlapped;
                    this->loc_.olap_length = overlap;
                    this->ctg_stats_.n_overlaps++;
                    this->ctg_stats_.length_overlap_tot += overlap;
                    this->ctg_stats_.length_olapd_loci_tot += fw_contig.length() + pe_contig.length() - overlap;
                    contig = move(fw_contig);
                    auto seq_itr = pe_contig.begin();
                    for(size_t i=0; i<size_t(overlap); ++i) {
                        assert(seq_itr != pe_contig.end());
                        ++seq_itr;
                    }
                    contig.append(seq_itr, pe_contig.end());
                } else {
                    if(detailed_output)
                        loc_.details_ss << "no_overlap\n";
                    this->loc_.ctg_status = LocData::separate;
                    contig = move(fw_contig);
                    for(size_t i=0; i<10; ++i)
                        contig.push_back(Nt4::n);
                    contig.append(pe_contig.begin(), pe_contig.end());
                }
            }
        }
    }
    return contig;
}

bool LocusProcessor::add_read_to_aln(
        CLocAlnSet& aln_loc,
        AlignRes& aln_res,
        SRead&& r,
        GappedAln* aligner,
        SuffixTree* stree
) {
    if (!this->align_reads_to_contig(stree, aligner, r.seq, aln_res))
        return false;

    Cigar cigar;
    parse_cigar(aln_res.cigar.c_str(), cigar);
    assert(!cigar.empty());
    simplify_cigar_to_MDI(cigar);
    cigar_extend_left(cigar, aln_res.subj_pos);
    assert(cigar_length_ref(cigar) <= aln_loc.ref().length());
    cigar_extend_right(cigar, aln_loc.ref().length() - cigar_length_ref(cigar));

    aln_loc.add(SAlnRead(move((Read&)r), move(cigar), r.sample));
    return true;
}

void LocusProcessor::phase_hets(
        vector<map<size_t,PhasedHet>>& phased_samples,
        vector<SiteCall>& calls,
        const CLocAlnSet& aln_loc,
        HaplotypeStats& hap_stats
) const {
    phased_samples.clear();
    phased_samples.resize(loc_.mpopi->samples().size());

    //
    // Apply a first MAC filter, if requested.
    //
    if (phasing_min_mac >= 2)
        for (SiteCall& call : calls)
            if (call.alleles().size() > 1)
                call.filter_mac(phasing_min_mac);

    //
    // Find SNPs.
    //
    vector<size_t> snp_cols;
    for (size_t i=0; i<aln_loc.ref().length(); ++i)
        if (calls[i].alleles().size() > 1)
            snp_cols.push_back(i);

    //
    // Check that the locus is polymorphic.
    //
    if (snp_cols.empty()) {
        // No SNPs, there is no phasing to do.
        ++hap_stats.n_badly_phased_samples[ {0, 0} ];
        return;
    }

    //
    // Initialize counters, streams.
    //
    size_t n_hets_needing_phasing = 0;
    size_t n_consistent_hets = 0;
    assert(hap_stats.per_sample_stats.size() == loc_.mpopi->samples().size());
    if (detailed_output)
        loc_.details_ss << "BEGIN phasing\n";
    stringstream o_hapgraph_ss;
    bool has_subgraphs = false;
    if (dbg_write_hapgraphs) {
        o_hapgraph_ss << "\n"
                      << "subgraph cluster_loc" << loc_.id << " {\n"
                      << "\tlabel=\"locus " << loc_.id << "\";\n"
                      << "\t# snp columns: ";
        join(snp_cols, ',', o_hapgraph_ss);
        o_hapgraph_ss << "\n";
    }

    //
    // Phase each sample.
    //
    vector<size_t> samples_failed;
    vector<size_t> samples_passed;
    for (size_t sample=0; sample<loc_.mpopi->samples().size(); ++sample) {
        if (phase_sample_hets(
                phased_samples[sample],
                calls, aln_loc, snp_cols, sample,
                hap_stats,
                n_hets_needing_phasing, n_consistent_hets,
                o_hapgraph_ss, has_subgraphs)
        ) {
            samples_passed.push_back(sample);
        } else {
            samples_failed.push_back(sample);
        }
    }

    //
    // Retry to phase failed samples after removing SNPs for which the minor
    // allele was only seen in hets and was never successfully phased.
    //
    if (!dbg_phasing_no_2ndpass) {
        vector<size_t> snp_cols_reduced;
        for(size_t col : snp_cols) {
            const vector<SampleCall>& c = calls[col].sample_calls();
            Counts<Nt2> counts;
            for (size_t sample : samples_passed) {
                switch(c[sample].call()) {
                case snp_type_hom:
                    counts.increment(c[sample].nt0(), 2);
                    break;
                case snp_type_het:
                    if (phased_samples[sample].count(col)) {
                        // This het was phased.
                        counts.increment(c[sample].nt0());
                        counts.increment(c[sample].nt1());
                    }
                    break;
                default:
                    break;
                }
            }
            size_t mac = counts.sorted()[1].first;
            if (mac >= size_t(phasing_min_mac))
                // Keep it.
                snp_cols_reduced.push_back(col);
        }
        if (snp_cols_reduced.size() < snp_cols.size())
            for (size_t sample : samples_failed)
                if (phase_sample_hets(
                        phased_samples[sample],
                        calls, aln_loc, snp_cols_reduced, sample,
                        hap_stats,
                        n_hets_needing_phasing, n_consistent_hets,
                        o_hapgraph_ss, has_subgraphs, false))
                    ++hap_stats.per_sample_stats[sample].n_phased_2ndpass;
    }

    //
    // Update the SiteCalls according to the haplotype calls, and reapply
    // the MAC filter on the updated data.
    //
    for (size_t col : snp_cols) {
        SiteCall& c = calls[col];
        for (size_t sample=0; sample<loc_.mpopi->samples().size(); ++sample)
            if (c.sample_calls()[sample].call() == snp_type_het
                    && !phased_samples[sample].count(col))
                c.discard_sample(sample);
        c.filter_mac(phasing_min_mac);
    }

    //
    // Clean up.
    //
    ++hap_stats.n_badly_phased_samples[ {n_consistent_hets, n_hets_needing_phasing} ];
    if (dbg_write_hapgraphs && has_subgraphs) {
        o_hapgraph_ss << "}\n";
        #pragma omp critical
        o_hapgraphs_f << o_hapgraph_ss.rdbuf();
    }
    if (detailed_output)
        loc_.details_ss << "END phasing\n";
}

bool LocusProcessor::phase_sample_hets(
        map<size_t,PhasedHet>& phased_sample,
        const vector<SiteCall>& calls,
        const CLocAlnSet& aln_loc,
        const vector<size_t>& snp_cols,
        size_t sample,
        HaplotypeStats& hap_stats,
        size_t& n_hets_needing_phasing, size_t& n_consistent_hets,
        ostream& o_hapgraph_ss, bool& has_subgraphs,
        bool first_pass
) const {
    phased_sample.clear();
    if (aln_loc.sample_reads(sample).empty())
        return true;
    if (first_pass)
        ++hap_stats.per_sample_stats[sample].n_diploid_loci;

    //
    // Find heterozygous positions.
    //
    vector<size_t> het_snps;
    for(size_t snp_i=0; snp_i<snp_cols.size(); ++snp_i)
        if (calls[snp_cols[snp_i]].sample_calls()[sample].call() == snp_type_het)
            het_snps.push_back(snp_i);
    if (!first_pass) {
        --hap_stats.per_sample_stats[sample].n_hets_2snps;
        --n_hets_needing_phasing;
    }
    if (het_snps.size() == 0) {
        // Sample is homozygous.
        return true;
    } else if (het_snps.size() == 1) {
        // One heterozygous SNP; sample has trivial 1nt-long haplotypes.
        size_t col = snp_cols[het_snps[0]];
        const SampleCall& c = calls[col].sample_calls()[sample];
        phased_sample.insert({col, {col, c.nt0(), c.nt1()}});
        return true;
    }
    ++hap_stats.per_sample_stats[sample].n_hets_2snps;
    ++n_hets_needing_phasing;

    vector<const SampleCall*> sample_het_calls; //**********
    for (size_t het_i=0; het_i<het_snps.size(); ++het_i)
        sample_het_calls.push_back(&calls[snp_cols[het_snps[het_i]]].sample_calls()[sample]);

    //
    // Iterate over reads, record as pairwise cooccurrences of  alleles.
    //
    SnpAlleleCooccurrenceCounter cooccurrences (snp_cols.size());
    count_pairwise_cooccurrences(cooccurrences, aln_loc,
                                    sample, snp_cols, het_snps, sample_het_calls);
    if (dbg_write_hapgraphs) {
        has_subgraphs = true;
        write_sample_hapgraph(o_hapgraph_ss, sample, het_snps, snp_cols,
                                sample_het_calls, cooccurrences);
    }

    //
    // Assemble phase sets.
    //
    vector<PhaseSet> phase_sets;
    bool phased = false;
    for (size_t min_n_cooccurrences = phasing_cooccurrences_thr_range.first;
            min_n_cooccurrences <= phasing_cooccurrences_thr_range.second;
            ++min_n_cooccurrences
    ) {
        if (assemble_phase_sets(
                phase_sets, het_snps, sample_het_calls,
                cooccurrences, min_n_cooccurrences
        )) {
            // Phasing succeeded.
            phased = true;
            if (detailed_output)
                loc_.details_ss << "phasing_ok\t"
                    << loc_.mpopi->samples()[sample].name
                    << "\tn_hets=" << het_snps.size()
                    << "\tcooc_thr=" << min_n_cooccurrences << '\n';
            break;
        }

        // Phasing failed. find alleles that are not connected to any
        // others, discard the het calls they are part of, and retry.
        if (phasing_dont_prune_hets)
            continue;
        vector<array<bool,2>> alleles_with_edges;
        alleles_with_edges.resize(het_snps.size(), {{false, false}});
        for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
            array<Nt2,2> allelesi = sample_het_calls[het_i]->nts();
            for (size_t het_j=het_i+1; het_j<het_snps.size(); ++het_j) {
                array<Nt2,2> allelesj = sample_het_calls[het_j]->nts();
                for (size_t allelei=0; allelei<2; ++allelei) {
                    for (size_t allelej=0; allelej<2; ++allelej) {
                        if (cooccurrences.at(het_snps[het_i],
                                                allelesi[allelei],
                                                het_snps[het_j],
                                                allelesj[allelej])
                                >= min_n_cooccurrences
                        ) {
                            alleles_with_edges[het_i][allelei] = true;
                            alleles_with_edges[het_j][allelej] = true;
                        }
                    }
                }
            }
        }
        vector<size_t> het_snps_reduced;
        vector<const SampleCall*> sample_het_calls_reduced;
        for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
            if (alleles_with_edges[het_i] == array<bool,2>({{true, true}})) {
                // Keep this het.
                het_snps_reduced.push_back(het_snps[het_i]);
                sample_het_calls_reduced.push_back(sample_het_calls[het_i]);
            }
        }
        if (het_snps_reduced.size() == het_snps.size())
            // No hets were removed; the situation is unchanged and phasing will fail again.
            continue;
        assert(!het_snps_reduced.empty()); // Otherwise the phasing wouldn't have previously failed.
        if (het_snps_reduced.size() == 1
            || assemble_phase_sets(phase_sets,
                                    het_snps_reduced,
                                    sample_het_calls_reduced,
                                    cooccurrences,
                                    min_n_cooccurrences)
        ) {
            // Phasing succeeded or became trivial.
            // We just downsize `het_snps`. We can do this because `phased_samples`
            // doesn't keep track of all the hets; only those of the best phase set.
            // In addition, unphased HETs are eventually discarded (as of Feb. 2018,
            // v2.0Beta8; c.f. PopMap::populate_locus) so we don't need to do anything
            // more to blacklist those genotypes.
            phased = true;
            if (detailed_output)
                loc_.details_ss << "phasing_ok_rmhets\t"
                    << loc_.mpopi->samples()[sample].name
                    << "\tn_hets=" << het_snps.size()
                    << "\tcooc_thr=" << min_n_cooccurrences
                    << "\tn_kept_hets=" << het_snps_reduced.size() << '\n';
            het_snps = het_snps_reduced;
            sample_het_calls = sample_het_calls_reduced;
            break;
        }
    }
    if (!phased) {
        if (detailed_output)
            loc_.details_ss << "phasing_failed\t"
                << loc_.mpopi->samples()[sample].name
                << "\tn_hets=" << het_snps.size() << '\n';
        return false;
    }
    ++hap_stats.per_sample_stats[sample].n_phased; // (if !first_pass, this wasn't reach before; no guard.)
    ++n_consistent_hets;
    if (het_snps.size() == 1) {
        // After pruning hets, sample has trivial 1nt-long haplotypes.
        // TODO Record that we've just removed the problem...
        size_t col = snp_cols[het_snps[0]];
        const SampleCall& c = *sample_het_calls[0];
        phased_sample.insert({col, {col, c.nt0(), c.nt1()}});
        return true;
    }

    //
    // Record phase sets.
    //
    for (const PhaseSet& ps : phase_sets) {
        assert(ps.size() == het_snps.size());
        size_t phase_set_id = SIZE_MAX;
        for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
            if (ps.het(het_i).left_allele == Nt4::n)
                continue;

            size_t col = snp_cols[het_snps[het_i]];
            if (phase_set_id == SIZE_MAX)
                phase_set_id = col; // Recommended value, c.f. VCF specification.
            PhasedHet ph = ps.het(het_i);
            ph.phase_set = phase_set_id;
            phased_sample[col] = ph;
        }
    }

    // Record singleton nodes (that are implicit in our representation).
    for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
        size_t col = snp_cols[het_snps[het_i]];
        if (!phased_sample.count(col)) {
            array<Nt2,2> alleles = sample_het_calls[het_i]->nts();
            phased_sample[col] = PhasedHet({col, alleles[0], alleles[1]});
        }
    }
    assert(!phased_sample.empty());

    //
    // Remove all phase sets but the largest one.
    //
    map<size_t,size_t> phase_set_sizes;
    for (auto& phasedhet : phased_sample)
        ++phase_set_sizes[phasedhet.second.phase_set];
    auto best_ps = phase_set_sizes.begin();
    for (auto ps=++phase_set_sizes.begin(); ps!=phase_set_sizes.end(); ++ps)
        if (ps->second > best_ps->second)
            best_ps = ps;
    for (auto phasedhet=phased_sample.begin(); phasedhet!=phased_sample.end();) {
        if (phasedhet->second.phase_set == best_ps->first)
            ++phasedhet;
        else
            phased_sample.erase(phasedhet++);
    }
    return true;
}

void LocusProcessor::count_pairwise_cooccurrences(
        SnpAlleleCooccurrenceCounter& cooccurrences,
        const CLocAlnSet& aln_loc,
        size_t sample,
        const vector<size_t>& snp_cols,
        const vector<size_t>& het_snps,
        const vector<const SampleCall*>& sample_het_calls
        ) const {

    cooccurrences.clear();
    vector<Nt4> read_hap (het_snps.size());
    for (size_t read_i : aln_loc.sample_reads(sample)) {
        // Build the haplotype.
        size_t curr_col = 0;
        Alignment::iterator read_itr (aln_loc.reads()[read_i].aln());
        for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
            size_t col = snp_cols[het_snps[het_i]];
            read_itr += col - curr_col;
            curr_col = col;
            Nt4 nt = *read_itr;
            if (nt == Nt4::n) {
                read_hap[het_i] = Nt4::n;
            } else {
                const SampleCall& c = *sample_het_calls[het_i];
                Nt2 nt2 = Nt2(nt);
                if (nt2 == c.nt0() || nt2 == c.nt1())
                    read_hap[het_i] = nt;
                else
                    read_hap[het_i] = Nt4::n;
            }
        }

        // Record the pairwise cooccurrences.
        for (size_t i=0; i<het_snps.size(); ++i) {
            Nt4 nti = read_hap[i];
            if (nti == Nt4::n)
                continue;
            for (size_t j=i+1; j<het_snps.size(); ++j) {
                Nt4 ntj = read_hap[j];
                if (ntj == Nt4::n)
                    continue;
                ++cooccurrences.at(het_snps[i], Nt2(nti), het_snps[j], Nt2(ntj));
            }
        }
    }
}

bool LocusProcessor::assemble_phase_sets(
        vector<PhaseSet>& phase_sets,
        const vector<size_t>& het_snps,
        const vector<const SampleCall*>& sample_het_calls,
        const SnpAlleleCooccurrenceCounter& cooccurrences,
        const size_t min_n_cooccurrences
) const {
    // xxx 2018 Feb. @Nick: Actually all of this would be more simple and efficient
    // if we used a single vector of PhasedHet's, with each het starting in its own
    // phase set.
    phase_sets.clear();

    // We keep track of which phase set each het is currently part of.
    vector<size_t> allele_to_ps (het_snps.size(), SIZE_MAX);

    assert(het_snps.size() == sample_het_calls.size());
    for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
        size_t snp_i = het_snps[het_i];
        size_t& ps_i = allele_to_ps[het_i]; //(by reference; always up to date)
        array<Nt2,2> allelesi = sample_het_calls[het_i]->nts();
        for (size_t het_j=het_i+1; het_j<het_snps.size(); ++het_j) {
            size_t snp_j = het_snps[het_j];
            size_t& ps_j = allele_to_ps[het_j];
            array<Nt2,2> allelesj = sample_het_calls[het_j]->nts();

            for (Nt2 nti : allelesi) {
                for (Nt2 ntj : allelesj) {
                    size_t n = cooccurrences.at(snp_i, nti, snp_j, ntj);
                    if (n < min_n_cooccurrences)
                        // Inexistent (n==0) or low-coverage edge.
                        continue;

                    if (ps_i == SIZE_MAX && ps_j == SIZE_MAX) {
                        // Both nodes are singletons. Start a new phase set.
                        ps_i = phase_sets.size();
                        ps_j = phase_sets.size();
                        phase_sets.push_back(PhaseSet(het_snps.size()));
                        phase_sets.back().add_het(het_i, allelesi);
                        phase_sets.back().add_het(het_j, allelesj, ntj, het_i, nti);
                    } else if (ps_i == SIZE_MAX) {
                        // Add `het_i` to `ps_j`.
                        phase_sets[ps_j].add_het(het_i, allelesi, nti, het_j, ntj);
                        ps_i = ps_j;
                    } else if (ps_j == SIZE_MAX) {
                        // Add `het_j` to `ps_i`.
                        phase_sets[ps_i].add_het(het_j, allelesj, ntj, het_i, nti);
                        ps_j = ps_i;
                    } else if (ps_i != ps_j) {
                        // Merge `ps_j` into `ps_i`.
                        phase_sets[ps_i].merge_with(phase_sets[ps_j], het_i, nti, het_j, ntj);
                        phase_sets[ps_j].clear();
                        size_t merged_ps = ps_j;
                        for (size_t& allele_ps : allele_to_ps)
                            if (allele_ps == merged_ps)
                                allele_ps = ps_i;
                    } else {
                        assert(ps_i == ps_j);
                        // Check that the edge is consistent.
                        if (!phase_sets[ps_i].is_edge_consistent(het_i, nti, het_j, ntj))
                            return false;
                    }
                }
            }
        }
    }

    // Purge empty phase sets.
    phase_sets.erase(std::remove_if(
            phase_sets.begin(), phase_sets.end(),
            [] (const PhaseSet& ps) {return ps.empty();}
            ),phase_sets.end());
    return true;
}

void PhaseSet::add_het(size_t het_i, array<Nt2,2> alleles) {
    assert(phased_hets_.size() > het_i);
    assert(phased_hets_ == PhaseSet(phased_hets_.size()).phased_hets_); // Uninitialized.

    phased_hets_[het_i].left_allele  = Nt4(alleles[0]);
    phased_hets_[het_i].right_allele = Nt4(alleles[1]);
}

void PhaseSet::add_het(size_t het_i, array<Nt2,2> alleles, Nt2 nt_i, size_t het_j, Nt2 nt_j) {
    assert(phased_hets_.size() > std::max(het_i, het_j));
    assert(phased_hets_[het_i].left_allele == Nt4::n); // `i` node shouldn't exist yet.
    assert(phased_hets_[het_j].left_allele != Nt4::n); //

    bool crossed_edge =
            (nt_j == phased_hets_[het_j].left_allele)
            ^ (nt_i == alleles[0]);

    if (crossed_edge) {
        // Flip the alleles so that the all the edges within the phase set are parallel.
        phased_hets_[het_i].left_allele  = Nt4(alleles[1]);
        phased_hets_[het_i].right_allele = Nt4(alleles[0]);
    } else {
        phased_hets_[het_i].left_allele  = Nt4(alleles[0]);
        phased_hets_[het_i].right_allele = Nt4(alleles[1]);
    }
}

void PhaseSet::merge_with(const PhaseSet& other, size_t het_i, Nt2 nt_i, size_t het_j, Nt2 nt_j) {
    assert(phased_hets_.size() > std::max(het_i, het_j));
    assert(phased_hets_.size() == other.phased_hets_.size());
    assert(phased_hets_[het_i].left_allele != Nt4::n);
    assert(other.phased_hets_[het_j].left_allele != Nt4::n);

    bool crossed_edge =
            (nt_i == phased_hets_[het_i].left_allele)
            ^ (nt_j == other.phased_hets_[het_j].left_allele);

    for (size_t i=0; i<phased_hets_.size(); ++i) {
        if (other.phased_hets_[i].left_allele != Nt4::n) {
            assert(other.phased_hets_[i].right_allele != Nt4::n); // Both alleles should be set/not set together.
            assert(phased_hets_[i].left_allele == Nt4::n); // Phase sets should be non-overlapping.
            if (crossed_edge) {
                phased_hets_[i].left_allele  = other.phased_hets_[i].right_allele;
                phased_hets_[i].right_allele = other.phased_hets_[i].left_allele;
            } else {
                phased_hets_[i].left_allele  = other.phased_hets_[i].left_allele;
                phased_hets_[i].right_allele = other.phased_hets_[i].right_allele;
            }
        }
    }
}

bool PhaseSet::is_edge_consistent(size_t het_i, Nt2 nt_i, size_t het_j, Nt2 nt_j) const {
    assert(phased_hets_[het_i].left_allele != Nt4::n);
    assert(phased_hets_[het_j].left_allele != Nt4::n);

    bool crossed_edge =
            (nt_i == phased_hets_[het_i].left_allele)
            ^ (nt_j == phased_hets_[het_j].left_allele);

    // We build phase sets so that all inside edges are parallel
    // (by flipping the alleles when necessary), so the proposed
    // edge is inconsistent if it is crossed.
    return !crossed_edge;
}

void LocusProcessor::write_one_locus (
        const CLocAlnSet& aln_loc,
        const vector<SiteCounts>& depths,
        const vector<SiteCall>& calls,
        const vector<map<size_t,PhasedHet>>& phase_data
) {
    char loc_id[16];
    sprintf(loc_id, "%d", loc_.id);

    const DNASeq4& ref = aln_loc.ref();
    const MetaPopInfo& mpopi = aln_loc.mpopi();

    //
    // Vcf output.
    //
    timers_.building_vcf.restart();
    assert(depths.size() == ref.length());
    assert(calls.size() == ref.length());
    vector<size_t> sample_sites_w_data (mpopi.samples().size(), 0);
    stringstream vcf_records;
    VcfRecord rec;
    for (size_t i=0; i<ref.length(); ++i) {
        const SiteCounts& sitedepths = depths[i];
        const SiteCall& sitecall = calls[i];
        if (sitecall.alleles().empty())
            // No useful data at this site.
            continue;

        if (!ref[i].is_acgt())
            continue;

        // Determine which alleles exist, and their order.
        // (n.b. As of Apr4,2017 the ref allele might not be the most frequent one.)
        vector<Nt2> vcf_alleles;
        map<Nt2, size_t> vcf_allele_indexes;
        {
            vcf_alleles.push_back(ref[i]);
            vcf_allele_indexes.insert({Nt2(ref[i]), 0});

            // Sort the alleles by frequency.
            vector<pair<double, Nt2>> sorted_alleles;
            for (auto& a : sitecall.alleles())
                sorted_alleles.push_back({a.second, a.first});
            sort(sorted_alleles.rbegin(), sorted_alleles.rend()); // (decreasing)

            // The reference allele has already been added to vcf_alleles; exclude it.
            for (auto iter=sorted_alleles.begin(); iter!=sorted_alleles.end(); ++iter) {
                if (iter->second == vcf_alleles[0]) {
                    sorted_alleles.erase(iter);
                    break;
                }
            }

            // Record the alternative alleles.
            for (auto& alt_allele : sorted_alleles) {
                vcf_allele_indexes.insert({alt_allele.second, vcf_alleles.size()});
                vcf_alleles.push_back(alt_allele.second);
            }
        }

        //
        // Create the VCF record.
        //

        // Chrom & pos.
        rec.clear();
        rec.append_chrom(loc_id);
        rec.append_pos(i);
        rec.append_id();

        // Alleles.
        for (Nt2 nt : vcf_alleles)
            rec.append_allele(nt);

        // SNP quality.
        rec.append_qual(sitecall.snp_qual());
        rec.append_filters();

        if(vcf_alleles.size() == 1) {
            // Fixed site.

            // Info/DP.
            rec.append_info(string("DP=") + to_string(sitedepths.tot.sum()));
            // Info/AD.
            Nt2 ref_nt = sitecall.alleles().begin()->first;
            rec.append_info(string("AD=") + to_string(sitedepths.tot[ref_nt]));
            if (dbg_write_nt_depths) {
                // Info/cnts.
                stringstream cnts;
                join(sitedepths.tot.arr(), ',', cnts);
                rec.append_info(string("cnts=") + cnts.str());
            }
            // Format.
            rec.append_format("DP");
            // Genotypes.
            for (size_t sample=0; sample<mpopi.samples().size(); ++sample) {
                size_t dp = sitedepths.samples[sample].sum();
                if (dp == 0) {
                    rec.append_sample(".");
                    continue;
                }
                ++sample_sites_w_data[sample];
                rec.append_sample(to_string(dp));
            }

        } else {
            // Polymorphic site.

            // Info/DP.
            rec.append_info(string("DP=") + to_string(sitedepths.tot.sum()));
            // Info/AD.
            vector<size_t> ad;
            for (Nt2 nt : vcf_alleles)
                ad.push_back(sitedepths.tot[nt]);
            stringstream ss;
            join(ad, ',', ss);
            rec.append_info(string("AD=") + ss.str());
            // Info/AF.
            vector<double> alt_freqs;
            for (auto nt=++vcf_alleles.begin(); nt!=vcf_alleles.end(); ++nt) // rem. always >1 alleles.
                alt_freqs.push_back(sitecall.alleles().at(*nt));
            rec.append_info(VcfRecord::util::fmt_info_af(alt_freqs));
            if (dbg_write_nt_depths) {
                // Info/cnts.
                stringstream cnts;
                join(sitedepths.tot.arr(), ',', cnts);
                rec.append_info(string("cnts=") + cnts.str());
            }

            // Format.
            rec.append_format("GT");
            if (!dbg_no_haplotypes) {
                rec.append_format("PS"); // Phase set.
                rec.append_format("FT"); // Filter.
            }
            rec.append_format("GQ");
            rec.append_format("DP");
            rec.append_format("AD");
            rec.append_format("GL");

            // Genotypes.
            for (size_t sample : mpopi.sample_indexes_orig_order()) {
                const Counts<Nt2>& sdepths = sitedepths.samples[sample];
                const SampleCall& scall = sitecall.sample_calls()[sample];
                if (sdepths.sum() == 0) {
                    // No data for this sample.
                    rec.append_sample("./.");
                    continue;
                }
                ++sample_sites_w_data[sample];

                stringstream genotype;
                // GT field.
                array<size_t,2> gt;
                switch (scall.call()) {
                case snp_type_hom:
                    gt[0] = vcf_allele_indexes.at(scall.nt0());
                    genotype << gt[0] << '/' << gt[0];
                    if (!dbg_no_haplotypes)
                        genotype << ":.:.";
                    genotype << ':' << scall.gq();
                    break;
                case snp_type_het:
                    if (!dbg_no_haplotypes) {
                        assert(!phase_data.empty());
                        const PhasedHet& p = phase_data[sample].at(i);
                        genotype << vcf_allele_indexes.at(p.left_allele)
                                    << '|'
                                    << vcf_allele_indexes.at(p.right_allele)
                                    << ':' << (p.phase_set + 1)
                                    << ':' << '.';
                    } else {
                        //dbg_no_haplotypes
                        gt[0] = vcf_allele_indexes.at(scall.nt0());
                        gt[1] = vcf_allele_indexes.at(scall.nt1());
                        std::sort(gt.begin(), gt.end()); // (Prevents '1/0'.)
                        genotype << gt[0] << '/' << gt[1];
                    }
                    genotype << ':' << scall.gq();
                    break;
                default:
                    genotype << "./.";
                    if (!dbg_no_haplotypes) {
                        if (scall.call() == snp_type_discarded)
                            genotype << ":.:disc";
                        else
                            genotype << ":.:.";
                    }
                    genotype << ":.";
                    break;
                }
                // DP field.
                genotype << ':' << sdepths.sum();
                // AD field.
                vector<size_t> ad;
                ad.reserve(vcf_alleles.size());
                for (Nt2 nt : vcf_alleles)
                    ad.push_back(sdepths[nt]);
                genotype << ':';
                join(ad, ',', genotype);
                // GL field.
                genotype << ':' << VcfRecord::util::fmt_gt_gl(vcf_alleles, scall.lnls());
                // cnts field.
                if (dbg_write_nt_depths) {
                    genotype << ":";
                    join(sdepths.arr(), ',', genotype);
                }
                // Push it.
                rec.append_sample(genotype.str());
            }
        }
        assert(rec.count_samples() == o_vcf_f->header().samples().size());

        vcf_records << rec;
    }
    loc_.o_vcf = vcf_records.str();
    timers_.building_vcf.update();

    //
    // Fasta output.
    //
    timers_.building_fa.restart();

    // Determine the number of samples for this locus. Some samples may have
    // been discarded (as of May 26, 2017, this would be because their haplotypes
    // were inconsistent).
    set<size_t> samples_w_reads;
    for (const SAlnRead& r : aln_loc.reads())
        samples_w_reads.insert(r.sample);
    size_t n_remaining_samples = 0;
    for (size_t sample_n_sites : sample_sites_w_data)
        if (sample_n_sites > 0)
            ++n_remaining_samples;


    // Write the record.
    string& fa = loc_.o_fa;
    assert(fa.empty());
    char cstr_buf[32];
    fa += '>';
    fa += loc_id;
    // Genomic position.
    if (!aln_loc.pos().empty()) {
        const PhyLoc& p = aln_loc.pos();
        sprintf(cstr_buf, "%u", p.bp+1);
        fa += " pos=";
        fa += p.chr();
        fa += ':';
        fa += cstr_buf;
        fa += ':';
        fa += (p.strand == strand_plus ? '+' : '-');
    }
    // Number of samples
    sprintf(cstr_buf, "%zu", n_remaining_samples);
    fa += " NS=";
    fa += cstr_buf;
    if (n_remaining_samples != samples_w_reads.size()) {
        assert(n_remaining_samples < samples_w_reads.size());
        sprintf(cstr_buf, "%zu", samples_w_reads.size() - n_remaining_samples);
        fa += " n_discarded_samples=";
        fa += cstr_buf;
    }
    // Contig status.
    switch (loc_.ctg_status) {
    case LocData::overlapped:
        sprintf(cstr_buf, "%zu", loc_.olap_length);
        fa += " contig=overlapped:";
        fa += cstr_buf;
        break;
    case LocData::separate:
        fa += " contig=separate";
        break;
    default:
        break;
    }
    fa += '\n';
    fa += ref.str();
    fa += '\n';
    timers_.building_fa.update();
}

void LocusProcessor::write_sample_hapgraph(
        ostream& os,
        size_t sample,
        const vector<size_t>& het_snps,
        const vector<size_t>& snp_cols,
        const vector<const SampleCall*>& sample_het_calls,
        const SnpAlleleCooccurrenceCounter& cooccurrences
        ) const {

    auto nodeid = [this,&sample](size_t col, Nt2 allele)
            {return string("l")+to_string(loc_.id)+"s"+to_string(sample)+"c"+to_string(col)+char(allele);};

    // Initialize the subgraph.
    os << "\tsubgraph cluster_sample" << sample << " {\n"
                  << "\t\tlabel=\""
                  << "i" << sample << " '" << loc_.mpopi->samples()[sample].name << "'\\n"
                  << "\";\n"
                  << "\t\tstyle=dashed;\n"
                  << "\t\t# heterozygous columns: ";
    vector<size_t> het_cols;
    for (size_t snp_i : het_snps)
        het_cols.push_back(snp_cols[snp_i]);
    join(het_cols, ',', os);
    os << "\n";

    // Write the node labels.
    for (size_t i=0; i<het_snps.size(); ++i) {
        array<Nt2,2> alleles = sample_het_calls[i]->nts();
        size_t col = snp_cols[het_snps[i]];
        for (Nt2 allele : alleles)
            os << "\t\t" << nodeid(col, allele)
                          << " [label=<"
                          << "<sup><font point-size=\"10\">" << col+1 << "</font></sup>" << allele
                          << ">];\n";
    }

    // Write the edges.
    for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
        array<Nt2,2> alleles_i = sample_het_calls[het_i]->nts();
        size_t snp_i = het_snps[het_i];
        for (size_t het_j=het_i+1; het_j<het_snps.size(); ++het_j) {
            array<Nt2,2> alleles_j = sample_het_calls[het_j]->nts();
            size_t snp_j = het_snps[het_j];
            for (Nt2 nti : alleles_i) {
                for (Nt2 ntj : alleles_j) {
                    size_t n = cooccurrences.at(snp_i, nti, snp_j, ntj);
                    if (n == 0)
                        continue;
                    os << "\t\t" << nodeid(snp_cols[snp_i],nti)
                                  << " -- " << nodeid(snp_cols[snp_j],ntj) << " [";
                    if (n==1)
                        os << "style=dotted";
                    else
                        os << "label=\"" << n << "\",penwidth=" << n;
                    os << "];\n";
                }
            }
        }
    }
    os << "\t}\n";
}

Timers& Timers::operator+= (const Timers& other) {
    reading += other.reading;
    processing += other.processing;
    writing_fa += other.writing_fa;
    writing_vcf += other.writing_vcf;
    writing_details += other.writing_details;
    writing_bams += other.writing_bams;
    // Within processing:
    processing_pre_alns += other.processing_pre_alns;
    rm_Ns += other.rm_Ns;
    assembling += other.assembling;
    init_alignments += other.init_alignments;
    aligning += other.aligning;
    building_bams += other.building_bams;
    merge_paired_reads += other.merge_paired_reads;
    processing_post_alns += other.processing_post_alns;
    rm_reads += other.rm_reads;
    genotyping += other.genotyping;
    haplotyping += other.haplotyping;
    building_vcf += other.building_vcf;
    building_fa += other.building_fa;
    cpt_consensus += other.cpt_consensus;
    counting_nts += other.counting_nts;
    return *this;
}

//
// Arguments
// ==========
//

const string help_string = string() +
        "gstacks " + VERSION  + "\n" +
        "\n"
        "De novo mode:\n"
        "  gstacks -P stacks_dir -M popmap\n"
        "\n"
        "  -P: input directory containg '*.matches.bam' files created by the\n"
        "      de novo Stacks pipeline, ustacks-cstacks-sstacks-tsv2bam\n"
        "\n"
        "Reference-based mode:\n"
        "  gstacks -I bam_dir -M popmap [-S suffix] -O out_dir\n"
        "  gstacks -B bam_file [-B ...] -O out_dir\n"
        "\n"
        "  -I: input directory containing BAM files\n"
        "  -S: with -I/-M, suffix to use to build BAM file names: by default this\n"
        "      is just '.bam', i.e. the program expects 'SAMPLE_NAME.bam'\n"
        "  -B: input BAM file(s)\n"
        "\n"
        "  The input BAM file(s) must be sorted by coordinate.\n"
        "  With -B, records must be assigned to samples using BAM \"reads groups\"\n"
        "  (gstacks uses the ID/identifier and SM/sample name fields). Read groups\n"
        "  must be consistent if repeated different files. Note that with -I, read\n"
        "  groups are unneeded and ignored.\n"
        "\n"
        "For both modes:\n"
        "  -M: path to a population map giving the list of samples\n"
        "  -O: output directory (default: none with -B; with -P same as the input\n"
        "      directory)\n"
        "  -t,--threads: number of threads to use (default: 1)\n"
        "\n"
        "SNP Model options:\n"
        "  --model: model to use to call variants and genotypes; one of\n"
        "           marukilow (default), marukihigh, or snp\n"
        "  --var-alpha: alpha threshold for discovering SNPs (default: 0.01 for marukilow)\n"
        "  --gt-alpha: alpha threshold for calling genotypes (default: 0.05)\n"
        "\n"
        "Paired-end options:\n"
        "  --rm-pcr-duplicates: remove all but one set ofread pairs of the same sample that \n"
        "                       have the same insert length (implies --rm-unpaired-reads)\n"
        "  --rm-unpaired-reads: discard unpaired reads\n"
        "  --ignore-pe-reads: ignore paired-end reads even if present in the input\n"
        "  --unpaired: ignore read pairing (only for paired-end GBS; treat READ2's as if they were READ1's)\n"
        "\n"
        "Advanced options:\n"
        "  (De novo mode)\n"
        "  --kmer-length: kmer length for the de Bruijn graph (default: 31, max. 31)\n"
        "  --max-debruijn-reads: maximum number of reads to use in the de Bruijn graph (default: 1000)\n"
        "  --min-kmer-cov: minimum coverage to consider a kmer (default: 2)\n"
        "  --write-alignments: save read alignments (heavy BAM files)\n"
        "\n"
        "  (Reference-based mode)\n"
        "  --min-mapq: minimum PHRED-scaled mapping quality to consider a read (default: 10)\n"
        "  --max-clipped: maximum soft-clipping level, in fraction of read length (default: 0.20)\n"
        "  --max-insert-len: maximum allowed sequencing insert length (default: 1000)\n"
        "\n"
        "  --details: write a heavier output\n"
        "  --phasing-cooccurrences-thr-range: range of edge coverage thresholds to\n"
        "        iterate over when building the graph of allele cooccurrences for\n"
        "        SNP phasing (default: 1,2)\n"
        "  --phasing-dont-prune-hets: don't try to ignore dubious heterozygote\n"
        "        genotypes during phasing\n"
        "\n"
#ifdef DEBUG
        "Debug options:\n"
        "  --dbg-no-overlaps: disable overlapping\n"
        "  --dbg-no-haps: disable phasing\n"
        "  --dbg-gfa: output a GFA file for each locus\n"
        "  --dbg-gfa-not-dag: output a GFA file for failed assemblies\n"
        "  --dbg-alns: output a file showing the contigs & alignments\n"
        "  --dbg-phasing-min-mac: minimum SNP MAC.\n"
        "  --dbg-phasing-no-2ndpass: don't try a second pass.\n"
        "  --dbg-hapgraphs: output a dot graph file showing phasing information\n"
        "  --dbg-depths: write detailed depth data in the output VCF\n"
        "  --dbg-log-stats-phasing: log detailed phasing statistics\n"
        "  --dbg-min-spl-reads: discard samples with less than this many reads (ref-based)\n"
        "  --dbg-min-loc-spls: discard loci with less than this many samples\n"
        "  --dbg-max-debruijn-reads\n"
        "\n"
#endif
        ;

string parse_command_line(int argc, char* argv[]) {
    auto bad_args = [](){
        cerr << help_string;
        exit(1);
    };
try {
    static const option long_options[] = {
        {"version",      no_argument,       NULL,  1000},
        {"help",         no_argument,       NULL,  'h'},
        {"quiet",        no_argument,       NULL,  'q'},
        {"stacks-dir",   required_argument, NULL,  'P'},
        {"in-dir",       required_argument, NULL,  'I'},
        {"in-bam",       required_argument, NULL,  'B'},
        {"suffix",       required_argument, NULL,  'S'},
        {"popmap",       required_argument, NULL,  'M'},
        {"out-dir",      required_argument, NULL,  'O'},
        {"unpaired",     no_argument,       NULL,  1007},
        {"threads",      required_argument, NULL,  't'},
        {"model",        required_argument, NULL,  1006},
        {"gt-alpha",     required_argument, NULL,  1005},
        {"var-alpha",    required_argument, NULL,  1008},
        {"kmer-length",  required_argument, NULL,  1001},
        {"max-debruijn-reads", required_argument, NULL, 1024},
        {"min-kmer-cov", required_argument, NULL,  1002},
        {"ignore-pe-reads", no_argument,    NULL,  1012},
        {"write-alignments", no_argument,   NULL,  1021},
        {"details",      no_argument,       NULL,  1013},
        {"min-mapq",     required_argument, NULL,  1014},
        {"max-clipped",  required_argument, NULL,  1015},
        {"max-insert-len", required_argument, NULL,  1016},
        {"rm-unpaired-reads", no_argument,  NULL,  1017},
        {"rm-pcr-duplicates", no_argument,  NULL,  1018},
        {"phasing-cooccurrences-thr-range", required_argument, NULL, 1019},
        {"phasing-dont-prune-hets", no_argument, NULL, 1020},
        //debug options
        {"min-kmer-freq", required_argument, NULL, 3021},
        {"dbg-phasing-min-mac", required_argument, NULL,  2018},
        {"dbg-phasing-no-2ndpass", no_argument, NULL, 2019},
        {"dbg-print-cloc-ids", no_argument, NULL,  2000},
        {"dbg-gfa",      no_argument,       NULL,  2003},
        {"dbg-gfa-not-dag",  no_argument,   NULL,  2015},
        {"dbg-alns",     no_argument,       NULL,  2004}, {"alns", no_argument, NULL, 2004},
        {"dbg-depths",   no_argument,       NULL,  2007},
        {"dbg-hapgraphs", no_argument,      NULL,  2010},
        {"dbg-no-overlaps", no_argument,    NULL,  2008},
        {"dbg-no-haps",  no_argument,       NULL,  2009},
        {"dbg-min-spl-reads", required_argument, NULL, 2014},
        {"dbg-min-loc-spls", required_argument, NULL, 2017},
        {0, 0, 0, 0}
    };

    string stacks_dir;
    string in_dir;
    string suffix = ".bam";

    double gt_alpha = 0.05;
    double var_alpha = 0.0;

    // bool pcr_dupl_measures = true;
    // auto pcr_duplicates_measures = [&](){
    //     pcr_dupl_measures = true;
    //     phasing_min_mac = 3;
    //     //TODO
    // };
    // auto no_pcr_duplicates_measures = [&](){
    //     pcr_dupl_measures = false;
    //     phasing_min_mac = 1;
    //     //TODO
    // };
    // pcr_duplicates_measures();

    int c;
    int long_options_i;
    std::cmatch tmp_cmatch;
    long tmp_long;
    while (true) {

        c = getopt_long(argc, argv, "hqP:I:B:S:s:M:O:W:t:m:", long_options, &long_options_i);

        if (c == -1)
            break;

        switch (c) {
        case 1000: //version
            cout << prog_name << " " << VERSION << "\n";
            exit(0);
            break;
        case 'h':
            cout << help_string;
            exit(0);
            break;
        case 'q':
            quiet = true;
            break;
        case 'P':
            stacks_dir = optarg;
            break;
        case 'I':
            in_dir = optarg;
            break;
        case 'B':
            in_bams.push_back(optarg);
            break;
        case 'S':
            suffix = optarg;
            break;
        case 'M':
            popmap_path = optarg;
            break;
        case 'O':
            out_dir = optarg;
            break;
        case 1007: //unpaired
            refbased_cfg.paired = false;
            break;
        case 1006: //model
            model_type = parse_model_type(optarg);
            break;
        case 1005: //gt-alpha
            gt_alpha = atof(optarg);
            break;
        case 1008: //var-alpha
            var_alpha = atof(optarg);
            break;
        case 't':
            num_threads = is_integer(optarg);
            if (num_threads < 0) {
                cerr << "Error: Illegal -t option value '" << optarg << "'.\n";
                bad_args();
            }
            break;
        case 1012: //ignore-pe-reads
            ignore_pe_reads = true;
            refbased_cfg.ign_pe_reads = true;
            break;
        case 1001://kmer-length
            km_length = atoi(optarg);
            if (km_length < 2 || km_length > 31) {
                cerr << "Error: Illegal -t option value '" << optarg << "' (valid range is 2-31).\n";
                bad_args();
            }
            break;
        case 1024: // max-debruijn-reads
            try {
                tmp_long = std::stol(optarg);
                if (tmp_long < 1)
                    throw exception();
            } catch(std::exception) {
                cerr << "Error: Illegal max-debruijn-reads option value '" << optarg << "'.\n";
                bad_args();
            }
            max_debruijn_reads = tmp_long;
            break;
        case 1002://min-kmer-cov
            min_km_count = atoi(optarg);
            break;
        case 1021://write-alignments
            bam_output = true;
            break;
        case 1013://details
            detailed_output = true;
            break;
        case 1014://min-mapq
            refbased_cfg.min_mapq = stoi(optarg);
            if (refbased_cfg.min_mapq > 255) {
                cerr << "Error: Illegal --min-mapq value '" << optarg << "'.\n";
                bad_args();
            }
            break;
        case 1015://max-clipped
            refbased_cfg.max_clipped = atof(optarg);
            if (refbased_cfg.max_clipped < 0.0 || refbased_cfg.max_clipped > 1.0) {
                cerr << "Error: Illegal --max-clipped value '" << optarg << "'.\n";
                bad_args();
            }
            break;
        case 1016://max-insert
            refbased_cfg.max_insert_refsize = stoi(optarg);
            break;
        case 1017://rm-unpaired-reads
            rm_unpaired_reads = true;
            break;
        case 1018://rm-pcr-duplicates
            rm_pcr_duplicates = true;
            rm_unpaired_reads = true;
            // no_pcr_duplicates_measures();
            break;
        // case 1022://pcr-duplicates-measures
        //     pcr_duplicates_measures();
        //     break;
        // case 1023://no-pcr-duplicates-measures
        //     no_pcr_duplicates_measures();
        //     break;
        case 1019: //phasing-cooccurrences-thr-range
            std::regex_search(optarg, tmp_cmatch, std::regex("^[0-9]+,[0-9]+$"));
            if (!tmp_cmatch.empty()) {
                phasing_cooccurrences_thr_range.first = atoi(optarg);
                phasing_cooccurrences_thr_range.second = atoi(strchr(optarg, ',') + 1);
            }
            if (tmp_cmatch.empty() || phasing_cooccurrences_thr_range.first > phasing_cooccurrences_thr_range.second) {
                cerr << "Error: Illegal --phasing-cooccurrences-thr-range value '"
                     << optarg << "'; expected 'INT_MIN,INT_MAX'.\n";
                bad_args();
            }
            break;
        case 1020: //phasing-dont-prune-hets
            phasing_dont_prune_hets = true;
            break;

        //
        // Debug options
        //
        case 3021://min-kmer-freq
            cerr << "WARNING: Ignored option --min-kmer-freq that does not exist anymore. See --max-debruijn-reads.\n";
            break;
        case 2018: //dbg-phasing-min-mac
            phasing_min_mac = is_integer(optarg);
            if (phasing_min_mac < 0) {
                cerr << "Error: Illegal --dbg-phasing-min-mac value '" << optarg << "'.\n";
                bad_args();
            }
            break;
        case 2019: //dbg-phasing-no-2ndpass
            dbg_phasing_no_2ndpass = true;
            break;
        case 2000://dbg-print-cloc-ids
            dbg_print_cloc_ids = true;
            break;
        case 2003://dbg-gfa
            dbg_write_gfa = true;
            break;
        case 2015://dbg-gfa-not-dag
            dbg_write_gfa_notdag = true;
            break;
        case 2004://dbg-alns
            dbg_write_alns = true;
            break;
        case 2010://dbg-hapgraphs
            dbg_write_hapgraphs = true;
            break;
        case 2007://dbg-depths
            dbg_write_nt_depths = true;
            break;
        case 2008://dbg-no-haps
            dbg_no_overlaps = true;
            break;
        case 2009://dbg-no-haps
            dbg_no_haplotypes = true;
            break;
        case 2014://dbg-min-spl-reads
            refbased_cfg.min_reads_per_sample = stoi(optarg);
            break;
        case 2017://dbg-min-loc-spls
            refbased_cfg.min_samples_per_locus = stoi(optarg);
            dbg_denovo_min_loc_samples = stoi(optarg);
            break;
        case '?':
            bad_args();
            break;
        default:
            DOES_NOT_HAPPEN;
            break;
        }
    }

    //
    // Check command consistency.
    //

    if (optind < argc) {
        cerr << "Error: Failed to parse command line: '" << argv[optind]
             << "' is seen as a positional argument. Expected no positional arguments.\n";
        bad_args();
    }

    size_t n_modes_given = !stacks_dir.empty() + !in_dir.empty() + !in_bams.empty();
    if (n_modes_given != 1) {
        cerr << "Error: Please specify exactly one of -P, -I or -B.\n";
        bad_args();
    }

    typedef GStacksInputT In;
    if (!stacks_dir.empty())
        input_type = popmap_path.empty()
            ? In::denovo_merger
            : In::denovo_popmap;
    else if (!in_dir.empty())
        input_type = In::refbased_popmap;
    else if (!in_bams.empty())
        input_type = In::refbased_list;
    else
        DOES_NOT_HAPPEN;

    if (input_type == In::refbased_popmap
        && popmap_path.empty()
    ) {
        cerr << "Error: Please specify a population map (-M).\n";
        bad_args();
    } else if (input_type == In::refbased_list
        && !popmap_path.empty()
    ) {
        cerr << "Error: Please specify -I/-M or -B, not both.\n";
        bad_args();
    }

    if (input_type == In::denovo_popmap || input_type == In::denovo_merger) {
        if (!refbased_cfg.paired) {
            cerr << "Error: --unpaired is for the reference-based mode.\n";
            bad_args();
        }
        if (suffix != ".bam") {
            cerr << "Error: --suffix is for the reference-based mode.\n";
            bad_args();
        }
    } else if (input_type == In::refbased_popmap || input_type == In::refbased_list) {
        if (!refbased_cfg.paired && rm_unpaired_reads) {
            cerr << "Error: --unpaired and --rm-unpaired-reads/--rm-pcr-duplicates are not compatible.\n";
            bad_args();
        }
        if (out_dir.empty()) {
            cerr << "Error: Please specify an output directory (-O).\n";
            bad_args();
        }
        if (min_km_count != 2 || max_debruijn_reads != 1000) {
            cerr << "Error: --min-kmer-count, --max-debruijn-reads are for the denovo mode.\n";
            bad_args();
        }
        if (bam_output) {
            cerr << "Error: --write-alignments is for the denovo mode.\n";
            bad_args();
        }
    } else {
        DOES_NOT_HAPPEN;
    }

    //
    // Process arguments.
    //

    for (string* dir : {&stacks_dir, &in_dir, &out_dir})
        if (!dir->empty() && dir->back() != '/')
            *dir += '/';

    if (var_alpha == 0.0) {
        if (model_type == marukilow) {
            var_alpha = 0.01;
        } else {
            cerr << "Error: No value was provided for --var-alpha"
                 << " (and there is no default for this model).\n";
            bad_args();
        }
    }
    switch (model_type) {
    case snp:        model.reset(new MultinomialModel(gt_alpha)); break;
    case marukihigh: model.reset(new MarukiHighModel(gt_alpha, var_alpha));  break;
    case marukilow:  model.reset(new MarukiLowModel(gt_alpha, var_alpha));   break;
    default:
        cerr << "Error: Model choice '" << to_string(model_type) << "' is not supported.\n";
        bad_args();
        break;
    }

    if (!popmap_path.empty()) {
        MetaPopInfo m;
        m.init_popmap(popmap_path);
        for (auto& s : m.samples())
            sample_names.push_back(s.name);
    }

    if (input_type == In::denovo_popmap) {
        if (out_dir.empty())
            out_dir = stacks_dir;
        for (const string& s : sample_names) {
            in_bams.push_back(stacks_dir + s + ".matches.bam");
            if (bam_output)
                out_bams.push_back(out_dir + s + ".alns.bam");
        }
    } else if (input_type == In::denovo_merger) {
        if (out_dir.empty())
            out_dir = stacks_dir;
        in_bams.push_back(stacks_dir + "catalog.bam");
        if (bam_output)
            out_bams.push_back(out_dir + "alignments.bam");
    } else if (input_type == In::refbased_popmap) {
        for (const string& sample : sample_names)
            in_bams.push_back(in_dir + sample + suffix);
    } else {
        assert(input_type == In::refbased_list);
    }

    check_or_mk_dir(out_dir);

    if (in_bams.empty())
        DOES_NOT_HAPPEN;

    //
    // Write the report.
    //
    stringstream os;
    os << "\n"
       << "Configuration for this run:\n";
    switch (input_type) {
    case GStacksInputT::denovo_popmap:
    case GStacksInputT::denovo_merger:
        os << "  Input mode: denovo\n";
        break;
    case GStacksInputT::refbased_popmap:
    case GStacksInputT::refbased_list:
        os << "  Input mode: reference-based";
        if (!refbased_cfg.paired)
            os << ", unpaired";
        os << "\n";
        break;
    default:
        DOES_NOT_HAPPEN;
        break;
    }
    if (!popmap_path.empty())
        os << "  Population map: '" << popmap_path << "'\n";
    os << "  Input files: " << in_bams.size() << ", e.g. '" << in_bams.front()<< "'\n"
       << "  Output to: '" << out_dir << "'\n"
       << "  Model: " << *model << "\n";
    if (ignore_pe_reads)
        os << "  Ignoring paired-end reads.\n";
    if (!refbased_cfg.paired)
        os << "  Ignoring pairing information.\n";
    if (km_length != 31)
        os << "  Kmer length: " << km_length << "\n";
    if (max_debruijn_reads != 1000)
        os << "  Maximum reads used for contig assembly: " << max_debruijn_reads << "\n";
    if (min_km_count != 2)
        os << "  Minimum kmer count: " << min_km_count << "\n";
    if (rm_unpaired_reads)
        os << "  Discarding unpaired reads.\n";
    if (rm_pcr_duplicates)
        os << "  Removing PCR duplicates.\n";
    // if (pcr_dupl_measures)
    //     os << "  PCR duplicates mitigation measures enabled.\n";

    return os.str();

} catch (std::invalid_argument&) {
    bad_args();
    exit(1);
}
}
