// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2017, Julian Catchen <jcatchen@illinois.edu>
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
#ifndef GSTACKS_H
#define GSTACKS_H

#include "constants.h"
#include "nucleotides.h"
#include "models.h"
#include "SuffixTree.h"
#include "GappedAln.h"
#include "BamI.h"
#include "locus_readers.h"

string parse_command_line(int argc, char* argv[]);

enum class GStacksInputT {unknown, denovo_popmap, denovo_merger, refbased_popmap, refbased_list};

//
// PhasedHet & PhaseSet
// ----------
//
struct PhasedHet {
    size_t phase_set; // N.B. The convention in VCF is to use the position of the first phased SNP for this.
    Nt4 left_allele;
    Nt4 right_allele;

    bool operator== (const PhasedHet& other) const
        {return phase_set == other.phase_set && left_allele == other.left_allele && right_allele == other.right_allele;}
};
class PhaseSet {
    vector<PhasedHet> phased_hets_;
public:
    PhaseSet() : phased_hets_{} {}
    void clear() {phased_hets_.clear();}
    bool empty() const {return phased_hets_.empty();}

    PhaseSet(size_t n_hets) : phased_hets_(n_hets, {SIZE_MAX, Nt4::n, Nt4::n}) {}
    size_t size() const {return phased_hets_.size();}
    const PhasedHet& het(size_t het_i) const {return phased_hets_[het_i];}

    // Add the first node to the phase set.
    void add_het(size_t het_i, array<Nt2,2> alleles);
    // Add a singleton node to the phase set, according to the given edge.
    void add_het(size_t het_i, array<Nt2,2> alleles, Nt2 nt_i, size_t het_j, Nt2 nt_j);
    // Merge with another phase set, according to the given edge.
    void merge_with(const PhaseSet& other, size_t het_i, Nt2 nt_i, size_t het_j, Nt2 nt_j);
    // Check that a new edge between two nodes within the phase set is consistent.
    bool is_edge_consistent(size_t het_i, Nt2 nt_i, size_t het_j, Nt2 nt_j) const;
};

//
// SnpAlleleCooccurrenceCounter
// ----------
// Convenience matrix class to record SNP allele co-occurences in
// reads/read pairs.
//
// The matrix is n_snps*n_snps but we only use the i<j half.
//
class SnpAlleleCooccurrenceCounter {
    size_t n_snps_;
    vector<array<array<size_t,4>,4>> cooccurrences_;
public:
    SnpAlleleCooccurrenceCounter(size_t n_snps) : n_snps_(n_snps), cooccurrences_(n_snps_*n_snps_) {};
    void clear();

    const size_t& at(size_t snp_i1, Nt2 snp1_allele, size_t snp_i2, Nt2 snp2_allele) const;
          size_t& at(size_t snp_i1, Nt2 snp1_allele, size_t snp_i2, Nt2 snp2_allele)
              {return (size_t&)((const SnpAlleleCooccurrenceCounter&)*this).at(snp_i1, snp1_allele, snp_i2, snp2_allele);}
};

//
// GenotypeStats
// ----------
// Statistics produced by LocusProcessor regarding genotyping.
//
class GenotypeStats {
public:
    size_t n_genotyped_loci;
    size_t n_sites_tot;
    double mean_n_sites_per_loc() const {return (double) n_sites_tot / n_genotyped_loci;}

    // Statistics for --rm-unpaired-reads, --rm-pcr-duplicates.
    struct PerSampleStats {
        size_t n_unpaired_reads;
        size_t n_read_pairs_pcr_dupl;
        size_t n_read_pairs_used;
        size_t n_loci_with_sample;
        size_t ns_cumsum;
        size_t ns_weighted_n_read_pairs_used;
    };
    vector<PerSampleStats> per_sample_stats;
    vector<size_t> pcr_clone_size_distrib;
    size_t n_unpaired_reads_rm() const {size_t n=0; for(auto& s: per_sample_stats) n+=s.n_unpaired_reads; return n;}
    size_t n_read_pairs_pcr_dupl() const {size_t n=0; for(auto& s: per_sample_stats) n+=s.n_read_pairs_pcr_dupl; return n;}
    size_t n_read_pairs_used() const {size_t n=0; for(auto& s: per_sample_stats) n+=s.n_read_pairs_used; return n;}

    GenotypeStats(size_t n_samples) : n_genotyped_loci(0), n_sites_tot(0), per_sample_stats(n_samples, PerSampleStats()) {}
    GenotypeStats& operator+= (const GenotypeStats& other);
};

//
// HaplotypeStats
// ----------
// Statistics produced by LocusProcessor regarding haplotyping.
//
class HaplotypeStats {
public:
    map<pair<size_t,size_t>,size_t> n_badly_phased_samples; // { {n_consistent_hets, n_hets_2snps} : count }

    struct PerSampleStats {
        size_t n_diploid_loci;
        size_t n_hets_2snps;
        size_t n_phased;
        size_t n_phased_2ndpass;
    };
    vector<PerSampleStats> per_sample_stats;

    HaplotypeStats(size_t n_samples) : per_sample_stats(n_samples, PerSampleStats()) {}
    HaplotypeStats& operator+= (const HaplotypeStats& other);

    size_t n_halplotyping_attempts() const;
    size_t n_consistent_hap_pairs() const;
};

//
// ContigStats
// ----------
// Statistics produced by `LocusProcessor::process()` when doing de novo assembly/alignment.
//
class ContigStats {
public:
    size_t n_nonempty_loci;
    size_t n_loci_w_pe_reads;
    size_t n_loci_almost_no_pe_reads;
    size_t n_loci_pe_graph_not_dag;
    size_t n_loci_pe_graph_fixed_to_dag;
    size_t length_ctg_tot;

    size_t n_aln_reads;
    size_t n_tot_reads;

    size_t n_overlaps;
    size_t length_overlap_tot;
    size_t length_olapd_loci_tot;

    OnlineMeanVar insert_length_olap_mv;

    size_t n_loci_no_pe_reads() const { return n_nonempty_loci - n_loci_w_pe_reads; }
    size_t n_loci_ctg() const { return n_loci_w_pe_reads - n_loci_almost_no_pe_reads - n_loci_pe_graph_not_dag; }
    double ctg_avg_length() const {return (double) length_ctg_tot / n_loci_ctg();}
    double mean_olap_length() const {return (double) length_overlap_tot / n_overlaps;}
    double mean_olapd_locus_length() const {return (double) length_olapd_loci_tot / n_overlaps;}

    ContigStats& operator+= (const ContigStats& other);
};

//
// LocData
// ----------
// Data regarding the locus LocusProcessor currently works on.
//
class LocData {
public:
    int id;
    PhyLoc pos;
    const MetaPopInfo* mpopi;

    enum ContigStatus {unknown, separate, overlapped};
    ContigStatus ctg_status = unknown;
    size_t olap_length;

    string o_vcf;
    string o_fa;
    string o_details;
    vector<vector<BamRecord>> o_bam;
    stringstream details_ss;

    void clear();
};

//
// Clocks
// ----------
// For Benchmarking. Thread specific.
//
struct Timers {
    Timer reading;
    Timer processing;
    Timer writing_fa;
    Timer writing_vcf;
    Timer writing_details;
    Timer writing_bams;

    Timer processing_pre_alns;
    Timer rm_Ns;
    Timer assembling;
    Timer init_alignments;
    Timer aligning;
    Timer building_bams;
    Timer merge_paired_reads;
    Timer processing_post_alns;
    Timer rm_reads;
    Timer genotyping;
    Timer haplotyping;
    Timer building_vcf;
    Timer building_fa;
    Timer cpt_consensus;
    Timer counting_nts;

    Timers& operator+= (const Timers& other);
};

//
// LocusProcessor
// ----------
// Functor for processing loci. Thread-specific.
//
class LocusProcessor {
public:
    LocusProcessor(size_t n_samples, const BamCLocReader* reader)
        : gt_stats_(n_samples), hap_stats_(n_samples), ctg_stats_(), loc_(), bam_cloc_reader_(reader) {}

    // Process a locus.
    void process(CLocReadSet& loc);
    void process(CLocAlnSet& aln_loc);

    // Access the output. Statistics & movable fasta/vcf per-locus text outputs.
    const GenotypeStats& gt_stats() const {return gt_stats_;}
    const HaplotypeStats& hap_stats() const {return hap_stats_;}
    const ContigStats& ctg_stats() const {return ctg_stats_;}
    Timers& timers() {return timers_;}

    string& vcf_out() {return loc_.o_vcf;}
    string& fasta_out() {return loc_.o_fa;}
    string& details_out() {return loc_.o_details;}
    vector<vector<BamRecord>>& bam_out() {return loc_.o_bam;}

private:
    GenotypeStats gt_stats_;
    HaplotypeStats hap_stats_;
    ContigStats ctg_stats_;
    Timers timers_;
    mutable LocData loc_;
    const BamCLocReader* bam_cloc_reader_;

    DNASeq4 assemble_locus_contig(
        const vector<SRead>& fw_reads,
        const vector<SRead>& pe_reads);
    bool add_read_to_aln(
            CLocAlnSet& aln_loc,
            AlignRes& aln_res,
            SRead&& read,
            GappedAln* aligner,
            SuffixTree* stree
            );

    bool align_reads_to_contig(SuffixTree *st, GappedAln *g_aln, DNASeq4 &query, AlignRes &aln_res) const;
    int suffix_tree_hits_to_dag(size_t query_len, vector<STAln> &alns, vector<STAln> &final_alns) const;
    size_t find_locus_overlap(SuffixTree *st, GappedAln *g_aln, const DNASeq4 &se_consensus, string &overlap_cigar) const;

    // Phase all samples.
    void phase_hets (
            vector<map<size_t,PhasedHet>>& phase_data,
            vector<SiteCall>& calls,
            const CLocAlnSet& aln_loc,
            HaplotypeStats& hap_stats
            ) const;

    // Phase (heterozygous) positions for one sample.
    bool phase_sample_hets(
            map<size_t,PhasedHet>& phased_sample,
            const vector<SiteCall>& calls,
            const CLocAlnSet& aln_loc,
            const vector<size_t>& snp_cols,
            size_t sample,
            HaplotypeStats& hap_stats,
            size_t& n_hets_needing_phasing,
            size_t& n_consistent_hets,
            ostream& o_hapgraph_ss,
            bool& has_subgraphs,
            bool first_pass=true) const;

    void count_pairwise_cooccurrences(
            SnpAlleleCooccurrenceCounter& cooccurrences,
            const CLocAlnSet& aln_loc,
            size_t sample,
            const vector<size_t>& snp_cols,
            const vector<size_t>& het_snps,
            const vector<const SampleCall*>& sample_het_calls
            ) const;

    // Assemble phase sets.
    // This is based on the graph of cooccurrences between the het alleles.
    // We build a graph in which the nodes are phased pairs of hets alleles
    // (`PhasedHet`s) and phase sets subgraphs of such nodes. The edges between
    // pairs of alleles are based on those in the graph of cooccurrences, but
    // in a phased context they appear as one of two types:
    // * parallel (+): left-allele-to-right-allele (left-to-left) or right-to-right.
    // * crossed (-): left-to-right, or right-to-left.
    // Allele linkage is transitive (through composition of +/- edges).
    // The assembly fails if there are inconsistent edges, i.e. edges that link
    // (transitively) one allele of one het SNP to both alleles of another het
    // SNP.
    // In practice, we start with an empty graph and iterate on edges, adding
    // nodes/PhasedHets when we first see them and checking that new edges are
    // consistent with the growing graph.
    bool assemble_phase_sets(
            vector<PhaseSet>& phase_sets,
            const vector<size_t>& het_snps,
            const vector<const SampleCall*>& sample_het_calls,
            const SnpAlleleCooccurrenceCounter& cooccurrences,
            const size_t min_n_cooccurrences
            ) const;

    // Create the fasta/vcf text outputs.
    void write_one_locus (
            const CLocAlnSet& aln_loc,
            const vector<SiteCounts>& depths,
            const vector<SiteCall>& calls,
            const vector<map<size_t,PhasedHet>>& phase_data // {col : phasedhet} maps, for all samples
    );

    // (debug) Write a sample's haplotype graph.
    void write_sample_hapgraph(
            ostream& os,
            size_t sample,
            const vector<size_t>& het_snps,
            const vector<size_t>& snp_cols,
            const vector<const SampleCall*>& sample_het_calls,
            const SnpAlleleCooccurrenceCounter& cooccurences
            ) const;

    //
    // (debug) Creates a CLocAlnSet from a CLocReadSet, using the true
    // reference but computing paired-end reads alignments. Same requirements as
    // `from_true_alignments()`.
    //
    void using_true_reference(CLocAlnSet& aln_loc, CLocReadSet&& loc);
};

#endif
