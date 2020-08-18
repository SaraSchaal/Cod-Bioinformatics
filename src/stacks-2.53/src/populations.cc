// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2012-2020, Julian Catchen <jcatchen@illinois.edu>
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
// populations -- generate population genetic statistics and output
// haplotypes in a population context.
//

#include <cctype>
#include <typeinfo>

#include "populations.h"
#include "export_formats.h"

using namespace std;

// Global variables to hold command-line options.
InputMode input_mode  = InputMode::stacks2;
int       num_threads =  1;
bool      quiet       = false;
string    in_path;
string    in_vcf_path;
string    out_path;
string    pmap_path;
string    bl_file;
string    wl_file;
string    bs_wl_file;
string    enz;
double    sigma             = 150000.0;
double    min_samples_overall   = 0.0;
double    min_samples_per_pop   = 0.0;
int       min_populations       = 1;
bool      filter_haplotype_wise = false;
int       batch_size        = 10000;
bool      calc_fstats       = false;
bool      calc_hwp          = false;
bool      bootstrap         = false;
bool      bootstrap_fst     = false;
bool      bootstrap_pifis   = false;
bool      bootstrap_phist   = false;
bool      bootstrap_div     = false;
bs_type   bootstrap_type    = bs_exact;
int       bootstrap_reps    = 100;
bool      bootstrap_wl      = false;
bool      write_single_snp  = false;
bool      write_random_snp  = false;
bool      merge_sites       = false;
bool      ordered_export    = false;
bool      smooth_fstats     = false;
bool      smooth_popstats   = false;
bool      loci_ordered      = false;
bool      mapping_cross     = false;
CrossT    mapcross_type     = CrossT::unk;
FormatT   mapcross_format   = FormatT::unk;
bool      log_fst_comp      = false;
bool      verbose           = false;
size_t    min_gt_depth      = 0;
double    merge_prune_lim   = 1.0;
double    minor_allele_freq = 0.0;
long      minor_allele_cnt  = 0;
double    max_obs_het       = 1.0;
double    p_value_cutoff    = 0.05;
corr_type fst_correction    = no_correction;
set<string> debug_flags;

string           out_prefix;
MetaPopInfo      mpopi;
set<int>         bootstraplist;

//
// Hold information about restriction enzymes
//
map<string, const char **> renz;
map<string, int>           renz_cnt;
map<string, int>           renz_len;
map<string, int>           renz_olap;

const int max_snp_dist = 500;

vector<Export *> exports;
template<typename E>
vector<Export*>::iterator find_export()
{
    return std::find_if(
            exports.begin(), exports.end(),
            [](Export* e){ return typeid(*e) == typeid(E); }
        );
}
template<typename E>
void add_export()
{
    // Check that the export isn't already present and add it.
    if (find_export<E>() == exports.end())
        exports.push_back(new E());
}

int main (int argc, char* argv[])
{
    unique_ptr<LogAlterator> logger;
try {

#ifndef HAVE_LIBZ
    cout << "Stacks was compiled without zlib, and will refuse to parse compressed files.\n";
#endif

    //
    // Initialize the globals that need it.
    //
    initialize_renz(renz, renz_cnt, renz_len);
    initialize_renz_olap(renz_olap);
    srandom(time(NULL));

    //
    // Parse the command line.
    //
    parse_command_line(argc, argv);

    //
    // Open and initialize the log file.
    //
    logger.reset(new LogAlterator(out_path + out_prefix, true, quiet, argc, argv));
    output_parameters(cout);

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    //
    // Read the population map file, if any.
    //
    if (!pmap_path.empty()) {
        cout << "Parsing population map...\n";
        mpopi.init_popmap(pmap_path);
        cout << "The population map contained " << mpopi.n_samples() << " samples, "
             << mpopi.pops().size() << " population(s), " << mpopi.groups().size() << " group(s).\n";
    } else {
        cout << "A population map was not specified, all samples will be read from '"
             << in_path << "' as a single popultaion.\n";
    }

    //
    // Locate and open input files, read VCF headers, parse population map, load black/white lists.
    //
    BatchLocusProcessor bloc(input_mode, batch_size, &mpopi);

    bloc.init(in_path, pmap_path);

    //
    // Initialize the genotypes processor If a mapping cross was specified.
    //
    MappingGenotypesProcessor *gproc = NULL;

    if (mapping_cross) {
        gproc = new MappingGenotypesProcessor(&mpopi, mapcross_type);
       
        // Check that the population map is compatible
        if (gproc->is_mapping_cross(cerr) == false)
            exit(1);

        //
        // Prepare outputs for the mapping cross data.
        //
        exports.push_back(new MarkersExport());

        switch (mapcross_format) {
        case FormatT::joinmap:
            exports.push_back(new JoinMapExport(gproc));
            break;
        case FormatT::rqtl:
            exports.push_back(new rQTLExport(gproc));
            break;
        case FormatT::onemap:
            if (mapcross_type == CrossT::cp)
                exports.push_back(new OneMapExport(gproc));
            else
                exports.push_back(new MapMakerExport(gproc));
            break;
        default:
            break;
        }
    }

    //
    // Report information on the structure of the populations specified.
    //
    mpopi.status(cout);

    if (mapping_cross == true && min_populations != 2) {
        cerr << "Notice: a mapping cross export has been specified, so the --min-populations parameter has been forced to '2'.\n\n";
        min_populations = 2;

    } else if (size_t(min_populations) > mpopi.pops().size()) {
        cerr << "Notice: --min-populations " << min_populations
             << " is larger than number of popualtions present, adjusting parameter to "
             << mpopi.pops().size() << "\n\n";
        min_populations = mpopi.pops().size();
    }

    //
    // Setup the default data exports.
    //
    exports.push_back(new RawHaplotypesExport());
    exports.push_back(new SumstatsExport());
    exports.push_back(new HapstatsExport());

    SnpDivergenceExport *sdiv_exp = NULL;
    HapDivergenceExport *hdiv_exp = NULL;
    if (calc_fstats) {
        sdiv_exp = new SnpDivergenceExport(logger->x);
        exports.push_back(sdiv_exp);
        hdiv_exp = new HapDivergenceExport();
        exports.push_back(hdiv_exp);
    }

    //
    // Setup the kernel smoothing apparatus.
    //
    LocusSmoothing *smooth = NULL;
    if (smooth_fstats || smooth_popstats)
        smooth = new LocusSmoothing(&mpopi, logger->x);

    //
    // Setup the divergence statistics calculator, if requested.
    //
    LocusDivergence *ldiv = NULL;
    if (calc_fstats)
        ldiv = new LocusDivergence(&mpopi);

    //
    // Open the export files and write any headers.
    //
    for (uint i = 0; i < exports.size(); i++) {
        exports[i]->open((const MetaPopInfo *) &mpopi);
        exports[i]->write_header();
    }

    //
    // Initialize the summary statistics object which will accumulate the metapopulation summary statistics
    // as the individual loci are read and processed.
    //
    SumStatsSummary sumstats(mpopi.pops().size());
    size_t n_sites = 0;
    size_t n_multiloci_sites = 0;

    const LocusFilter &filter = bloc.filter();
    cout << "\nProcessing data in batches:\n"
         << "  * load a batch of catalog loci and apply filters\n"
         << "  * compute SNP- and haplotype-wise per-population statistics\n";
    if (calc_hwp)
        cout << "  * compute SNP- and haplotype-wise deviation from HWE\n";
    if (calc_fstats)
        cout << "  * compute F-statistics\n";
    if (smooth_popstats)
        cout << "  * smooth per-population statistics\n";
    if (smooth_fstats)
        cout << "  * smooth F-statistics\n";
    cout << "  * write the above statistics in the output files\n"
         << "  * export the genotypes/haplotypes in specified format(s)\n"
         << "More details in '" << logger->distribs_path << "'.\n"
         << "Now processing...\n"
         << flush;

    Timer timer;
    logger->x << "\nBEGIN batch_progress\n";
    while (true) {
        //
        // Read the next set of loci to process.
        // - If data are denovo, load blim._batch_size loci.
        // - If data are reference aligned, load one chromosome.
        // - Filter the loci according to command line parameters (-r, -p, --maf, --write-single-snp, etc.)
        // - Sort the loci by basepair if they are ordered.
        //
        size_t loc_cnt = bloc.next_batch(logger->x);
        if (filter.batch_seen() == 0)
            break;

        if (loci_ordered) {
            cout << bloc.chr() << " " << flush;
            logger->x << bloc.chr();
        } else {
            cout << "Batch " << bloc.next_batch_number() - 1 << " " << flush;
            logger->x << "Batch " << bloc.next_batch_number() - 1;
        }
        logger->x << ": analyzed "
                  << filter.batch_total() << " loci; filtered "
                  << filter.batch_filtered() << " loci; "
                  << filter.batch_seen() << " loci seen.\n";

        if (loc_cnt == 0)
            break;
        assert(!bloc.loci().empty());

        sumstats.accumulate(bloc.loci());

        //
        // Calculate haplotype and gene diversity, Hardy-Weinberg proportions per locus per population.
        //
        bloc.hapstats(logger->x);

        //
        // Calculate and report the extent of overlap between different RAD loci.
        //
        if (loci_ordered) {
            size_t chr_n_sites;
            size_t chr_n_multiloci_sites;
            bloc.report_locus_overlap(chr_n_sites, chr_n_multiloci_sites, (verbose ? &logger->x : NULL));
            n_sites += chr_n_sites;
            n_multiloci_sites += chr_n_multiloci_sites;
            logger->x << "    " << chr_n_sites << " genomic sites, of which "
                      << chr_n_multiloci_sites <<  " were covered by multiple loci ("
                      << as_percentage(chr_n_multiloci_sites, chr_n_sites) << ").\n";
        }

        //
        // Calculate divergence statistics (Fst), if requested.
        //
        if (calc_fstats) {
            logger->x << "    Calculating F statistics...";
            ldiv->snp_divergence(bloc.loci());
            ldiv->haplotype_divergence_pairwise(bloc.loci());
            ldiv->haplotype_divergence(bloc.loci());
            logger->x << "done.\n";
        }

        //
        // Smooth population statistics across individual populations, and between populations.
        //
        if ( (smooth_fstats || smooth_popstats) && loci_ordered == false) {
            logger->x << "    Notice: Smoothing was requested (-k), but will not be performed as the loci are not ordered.\n";
        } else if (smooth_fstats || smooth_popstats) {
            logger->x << "    Generating kernel-smoothed population statistics...";
            if (smooth_popstats) {
                smooth->snpstats(bloc.loci(), logger->x);
                smooth->hapstats(bloc.loci(), logger->x);
            }
            if (smooth_fstats) {
                smooth->snp_divergence(bloc.loci(), ldiv->snp_values(), logger->x);
                smooth->hap_divergence(bloc.loci(), ldiv->haplotype_values(), ldiv->metapop_haplotype_values(), logger->x);
            }
            logger->x << "done.\n";
        }

        //
        // Process markers for genetic mapping cross.
        //
        if (mapping_cross)
            gproc->next_batch(bloc.loci());

        //
        // Export this subset of the loci.
        //
        for (uint i = 0; i < exports.size(); i++)
            exports[i]->write_batch(bloc.loci());

        if (calc_fstats) {
            sdiv_exp->write_batch_pairwise(bloc.loci(), ldiv->snp_values());
            hdiv_exp->write_batch_pairwise(bloc.loci(), ldiv->haplotype_values(), ldiv->metapop_haplotype_values());
            ldiv->clear(bloc.loci());
        }

        logger->x << flush;
        timer.update();
        #ifdef DEBUG
        cout << "(" << (size_t) timer.elapsed() << "s)\n" << flush;
        #else
        cout << "\n";
        #endif

    }
    logger->x << "END batch_progress\n";

    //
    // Report what we read from the input files.
    //
    bloc.summarize(cout);

    cout << "\n"
         << "Removed " << filter.filtered() << " loci that did not pass sample/population constraints from " << filter.seen() << " loci.\n"
         << "Kept " << filter.total() << " loci, composed of " << filter.total_sites() << " sites; "
         << filter.filtered_sites() << " of those sites were filtered, " << filter.variant_sites() << " variant sites remained.\n";
    if (loci_ordered)
        cout << "    " << n_sites << " genomic sites, of which "
             << n_multiloci_sites <<  " were covered by multiple loci ("
             << as_percentage(n_multiloci_sites, n_sites) << ").\n";
    if (filter.total_sites() == 0) {
        cerr << "Error: All data has been filtered out.\n";
        throw exception();
    }

    //
    // Do the final sumstats calculations and write the sumstats summary files.
    //
    sumstats.final_calculation();
    sumstats.write_results();

    if (calc_hwp) {
        cout << "Number of loci found to be significantly out of Hardy-Weinberg equilibrium (<" << p_value_cutoff << "):\n";
        for (uint j = 0; j < mpopi.pops().size(); j++)
            cout << "  "
                 << mpopi.pops()[j].name << ": "
                 << bloc._sig_hwe_dev[j] << "\n";
        cout << "(more detail in populations.sumstats.tsv and populations.hapstats.tsv)\n";
    }

    if (calc_fstats)
        ldiv->write_summary(out_path + out_prefix);

    //
    // Write out the distributions of catalog loci.
    //
    bloc.write_distributions(logger->x);

    if (smooth_fstats || smooth_popstats)
        delete smooth;
    if (calc_fstats)
        delete ldiv;

    //
    // Close the export files and do any required post processing.
    //
    for (uint i = 0; i < exports.size(); i++) {
        exports[i]->post_processing();
        exports[i]->close();
        delete exports[i];
    }

    if (mapping_cross)
        delete gproc;
    
    cout << "Populations is done.\n";
    return 0;

} catch (exception& e) {
    return stacks_handle_exceptions(e);
}
}

int
BatchLocusProcessor::init(string in_path, string pmap_path)
{
    //
    // Read the blacklist and whitelist to control which loci we load.
    //
    int cnt;
    if (bl_file.length() > 0) {
        cnt = this->_loc_filter.load_blacklist(bl_file);
        cout << "Loaded " << cnt << " blacklisted markers.\n";
    }
    if (wl_file.length() > 0) {
        cnt = this->_loc_filter.load_whitelist(wl_file);
        cout << "Loaded " << cnt << " whitelisted markers.\n";
        this->_user_supplied_whitelist = true;
    }

    if (this->_input_mode == InputMode::vcf)
        this->init_external_loci(in_vcf_path, pmap_path);
    else
        this->init_stacks_loci(in_path, pmap_path);

    //
    // Initialize our per-population haplotype counters.
    //
    this->_sig_hwe_dev = new size_t[this->_mpopi->pops().size()];
    memset(this->_sig_hwe_dev, 0, sizeof(size_t)*this->_mpopi->pops().size());

    return 0;
}

size_t
BatchLocusProcessor::next_batch(ostream &log_fh)
{
    this->_batch_num++;
    if (this->_input_mode == InputMode::vcf)
        return this->next_batch_external_loci(log_fh);
    else
        return this->next_batch_stacks_loci(log_fh);
}

void BatchLocusProcessor::batch_clear()
{
    for (LocBin* loc : this->_loci)
        delete loc;
    this->_loci.clear();
    this->_chr.clear();
    this->_loc_filter.batch_clear();
}

int
BatchLocusProcessor::init_stacks_loci(string in_path, string pmap_path)
{
    //
    // Open the files.
    //
    string catalog_fa_path  = in_path + "catalog.fa.gz";
    string catalog_vcf_path = in_path + "catalog.calls";

    this->_fasta_reader.open(catalog_fa_path);
    this->_cloc_reader.open(catalog_vcf_path);

    // Create the population map or check that all samples have data.
    if (pmap_path.empty()) {
        cout << "No population map specified, using all samples...\n";
        this->_mpopi->init_names(this->cloc_reader().header().samples());
    } else {
        size_t n_samples_before = this->_mpopi->n_samples();
        this->_mpopi->intersect_with(this->cloc_reader().header().samples());
        size_t n_rm_samples = n_samples_before - this->_mpopi->n_samples();

        if (n_rm_samples > 0) {
            cerr << "Warning: No genotype data exists for " << n_rm_samples
                 << " of the samples listed in the population map.\n";
            if (this->_mpopi->samples().empty()) {
                cerr << "Error: No more samples.\n";
                throw exception();
            }
        }
    }

    this->cloc_reader().set_sample_ids(*this->_mpopi);

    //
    // Initialize the locus filter after we have constructed the population map.
    //
    this->_loc_filter.init(this->pop_info());

    //
    // Create the {sample_vcf_index : sample_popmap_index} table.
    //
    const vector<string>& vcf_samples = this->_cloc_reader.header().samples();
    this->_samples_vcf_to_mpopi.assign(vcf_samples.size(), SIZE_MAX);
    for (size_t vcf_i=0; vcf_i<vcf_samples.size(); ++vcf_i) {
        size_t mpopi_i = this->_mpopi->get_sample_index(vcf_samples[vcf_i], false);
        if (mpopi_i != SIZE_MAX)
            this->_samples_vcf_to_mpopi[vcf_i] = mpopi_i;
    }

    return 0;
}

size_t
BatchLocusProcessor::next_batch_stacks_loci(ostream &log_fh)
{
    this->batch_clear();

    //
    // Check if we queued a LocBin object from the last round of reading.
    //
    if (this->_next_loc != NULL) {
        this->_loci.push_back(this->_next_loc);
        this->_next_loc = NULL;
        this->_loc_filter.locus_seen();
        this->_loc_filter.keep_locus(this->_loci.back());
    }

    Seq seq;
    seq.id = new char[id_len];
    vector<VcfRecord> records;
    while(this->_cloc_reader.read_one_locus(records)) {
        this->_loc_filter.locus_seen();

        //
        // Get the current locus ID and find the corresponding fasta record.
        // (Note: c-loci with very low coverage might be entirely missing from
        // the VCF; in this case ignore them.)
        //
        assert(!records.empty());
        int cloc_id = is_integer(records[0].chrom());
        assert(cloc_id >= 0);
        do {
            int rv = this->_fasta_reader.next_seq(seq);
            if (rv == 0) {
                cerr << "Error: catalog VCF and FASTA files are discordant, maybe trucated. rv: "
                     << rv << "; cloc_id: " << cloc_id << "\n";
                throw exception();
            }
        } while (is_integer(seq.id) != cloc_id);

        //
        // Check if this locus is white/blacklisted.
        //
        if (this->_loc_filter.whitelist_filter(cloc_id) ||
            this->_loc_filter.blacklist_filter(cloc_id))
            continue;

        //
        // Create and populate a new catalog locus & the associated genotypes.
        //
        LocBin* loc = new LocBin(this->_mpopi->n_samples());
        loc->cloc = new CSLocus();
        loc->d = new Datum *[this->_mpopi->n_samples()];
        for (size_t i=0; i<this->_mpopi->n_samples(); ++i)
            loc->d[i] = NULL;

        PopMap<CSLocus>::populate_internal(
            loc->cloc, loc->d,
            seq, records,
            this->_cloc_reader.header(), this->_mpopi, this->_samples_vcf_to_mpopi);

        //
        // Apply locus & SNP filters.
        //
        this->_dists.accumulate_pre_filtering(loc->sample_cnt, loc->cloc);
        this->_loc_filter.whitelist_snp_filter(*loc);
        if (this->_loc_filter.apply_filters_stacks(*loc, log_fh, *this->_mpopi)) {
            delete loc;
            continue;
        }
        assert(loc->s != NULL);

        //
        // Detect the end of batch.
        //
        loci_ordered = !loc->cloc->loc.empty();
        if (loci_ordered) {
            //
            // Ref-based.
            //
            if (this->_chr.empty()) {
                this->_chr = loc->cloc->loc.chr();
            } else if (this->_chr.compare(loc->cloc->loc.chr()) != 0) {
                if (!this->_loci.empty()) {
                    this->_next_loc = loc;
                    this->_loc_filter.locus_unsee();
                    break;
                } else {
                    //
                    // No loci were analyzed on the previous chromosome, keep processing
                    // the current chromosome without returning.
                    //
                    this->_chr = loc->cloc->loc.chr();
                }
            }
            this->_loci.push_back(loc);
            this->_loc_filter.keep_locus(loc);
        } else {
            //
            // De novo.
            //
            if (ordered_export) {
                cerr << "Error: Options --ordered-export/--smooth/--bootstrap"
                    << " are for reference-aligned data only.\n";
                throw exception();
            }
            loc->cloc->loc.set("un", this->_unordered_bp, strand_plus);
            this->_unordered_bp += loc->cloc->len;
            this->_loci.push_back(loc);
            this->_loc_filter.keep_locus(loc);
            if (this->_loci.size() == this->_batch_size)
                break;
        }
    }

    //
    // Record the post-filtering distribution of catalog loci for this batch.
    //
    this->_dists.accumulate(this->_loci);

    //
    // Sort the catalog loci, if possible.
    //
    if (loci_ordered)
        sort(this->_loci.begin(), this->_loci.end(),
             [] (const LocBin *a, const LocBin *b) -> bool {
                 return a->cloc->loc.bp < b->cloc->loc.bp;
             });

    return this->_loci.size();
}

int
BatchLocusProcessor::init_external_loci(string in_path, string pmap_path)
{
    //
    // Open the VCF file
    //
    cout << "Opening the VCF file...\n";
    this->_vcf_parser.open(in_path);

    if (this->_vcf_parser.header().samples().empty()) {
        cerr << "Error: No samples in VCF file '" << in_path << "'.\n";
        throw exception();
    }

    // Reconsider the MetaPopInfo in light of the VCF header.
    if (pmap_path.empty()) {
        cout << "No population map specified, creating one from the VCF header...\n";
        this->_mpopi->init_names(this->_vcf_parser.header().samples());

    } else {
        // Intersect the samples present in the population map and the VCF.
        size_t n_samples_before = this->_mpopi->n_samples();

        this->_mpopi->intersect_with(this->_vcf_parser.header().samples());

        size_t n_rm_samples = n_samples_before - this->_mpopi->n_samples();
        if (n_rm_samples > 0) {
            cerr << "Warning: Of the samples listed in the population map, "
                 << n_rm_samples << " could not be found in the VCF :";
            if (this->_mpopi->samples().empty()) {
                cerr << "Error: No more samples.\n";
                throw exception();
            }
        }
    }

    // Create arbitrary sample IDs.
    for (size_t i = 0; i < this->_mpopi->n_samples(); ++i)
        this->_mpopi->set_sample_id(i, i+1); //id=i+1

    //
    // Initialize the locus filter after we have constructed the population map.
    //
    this->_loc_filter.init(this->pop_info());

    //
    // Create the {sample_vcf_index : sample_popmap_index} table.
    //
    const vector<string>& vcf_samples = this->_vcf_parser.header().samples();
    this->_samples_vcf_to_mpopi.assign(vcf_samples.size(), SIZE_MAX);
    for (size_t vcf_i=0; vcf_i<vcf_samples.size(); ++vcf_i) {
        size_t mpopi_i = this->_mpopi->get_sample_index(vcf_samples[vcf_i], false);
        if (mpopi_i != SIZE_MAX)
            this->_samples_vcf_to_mpopi[vcf_i] = mpopi_i;
    }

    this->_total_ext_vcf = 0;

    return 0;
}

size_t
BatchLocusProcessor::next_batch_external_loci(ostream &log_fh)
{
    //
    // VCF mode
    //
    this->batch_clear();
    loci_ordered = true;

    //
    // Check if we queued a LocBin object from the last round of reading.
    //
    if (this->_next_loc != NULL) {
        this->_loci.push_back(this->_next_loc);
        this->_next_loc = NULL;
        this->_loc_filter.locus_seen();
        this->_loc_filter.keep_locus(this->_loci.back());
    }

    int cloc_id = (this->_loci.empty() ? 1 : this->_loci.back()->cloc->id + 1);
    VcfRecord rec;
    while (this->_vcf_parser.next_record(rec)) {
        this->_loc_filter.locus_seen();
        this->_total_ext_vcf++;

        // Check for a SNP.
        if (not rec.is_snp()) {
            this->_skipped_notsnp.push_back(this->_vcf_parser.line_number());
            continue;
        } else if (rec.count_alleles() != 2) {
            this->_skipped_notbinarysnp.push_back(this->_vcf_parser.line_number());
            continue;
        }

        // Check for a filtered-out SNP
        if (strncmp(rec.filters(), ".", 2) != 0 && strncmp(rec.filters(), "PASS", 5) != 0) {
            this->_skipped_filter.push_back(this->_vcf_parser.line_number());
            continue;
        }

        //
        // Create and populate a new catalog locus.
        //
        LocBin* loc = new LocBin(this->_mpopi->n_samples());
        loc->cloc = new CSLocus();
        loc->d = new Datum *[this->_mpopi->n_samples()];
        for (size_t i = 0; i < this->_mpopi->n_samples(); i++)
            loc->d[i] = NULL;
        if (!PopMap<CSLocus>::populate_external(
                loc->cloc, loc->d,
                cloc_id++, rec,
                this->_vcf_parser.header(), this->_mpopi, this->_samples_vcf_to_mpopi)
        ) {
            // Bad record; a warning has been printed.
            delete loc;
            continue;
        }
        assert(!loc->cloc->loc.empty());
        assert(loc->cloc->len == 1);
        assert(loc->cloc->snps.size() == 1);

        //
        // Apply filters.
        //
        this->_dists.accumulate_pre_filtering(loc->sample_cnt, loc->cloc);
        if (this->_loc_filter.apply_filters_external(*loc, log_fh, *this->_mpopi)) {
            delete loc;
            continue;
        }
        assert(loc->s != NULL);

        //
        // Detect the end of batch.
        //
        if (this->_chr.empty()) {
            this->_chr = loc->cloc->loc.chr();
        } else if (this->_chr.compare(loc->cloc->loc.chr()) != 0) {
            this->_next_loc = loc;
            this->_loc_filter.locus_unsee();
            break;
        }
        this->_loci.push_back(loc);
        this->_loc_filter.keep_locus(loc);
    }

    //
    // Record the post-filtering distribution of catalog loci for this batch.
    //
    this->_dists.accumulate(this->_loci);

    //
    // Sort the catalog loci, if possible.
    //
    sort(this->_loci.begin(), this->_loci.end(),
        [] (const LocBin *a, const LocBin *b) -> bool {
            return a->cloc->loc.bp < b->cloc->loc.bp;
        });

    return this->_loci.size();
}

int
BatchLocusProcessor::summarize(ostream &log_fh)
{
    if (this->_input_mode == InputMode::vcf) {
        log_fh << "Found " << this->_total_ext_vcf << " SNP records in file '" << in_vcf_path
             << "'. (Skipped " << this->_skipped_filter.size() << " already filtered-out SNPs and "
             << this->_skipped_notsnp.size() << " non-SNP records ; more with --verbose.)\n";
        if (verbose && not this->_skipped_notsnp.empty()) {
            log_fh << "The following VCF record lines were determined not to be SNPs and skipped :";
            for (vector<size_t>::const_iterator l = this->_skipped_notsnp.begin(); l != this->_skipped_notsnp.end(); ++l)
                log_fh << " " << *l;
            log_fh << "\n";
        }
        if (verbose && not this->_skipped_notbinarysnp.empty()) {
            log_fh << "The following VCF record lines were determined not to be binary SNPs and skipped :";
            for (vector<size_t>::const_iterator l = this->_skipped_notbinarysnp.begin(); l != this->_skipped_notbinarysnp.end(); ++l)
                log_fh << " " << *l;
            log_fh << "\n";
        }
    } else {
    }

    return 0;
}

int
BatchLocusProcessor::hapstats(ostream &log_fh)
{
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (uint i = 0; i < this->_loci.size(); i++) {
            LocBin *loc = this->_loci[i];

            loc->s->calc_hapstats(loc->cloc, (const Datum **) loc->d, *this->_mpopi);
        }
    }

    if (calc_hwp) {
        const LocStat *l;
        const vector<Pop> &pops = this->_mpopi->pops();

        for (uint i = 0; i < this->_loci.size(); i++) {
            for (uint j = 0; j < pops.size(); j++) {
                l = this->_loci[i]->s->hapstats_per_pop(j);

                if (l == NULL)
                    continue;

                this->_sig_hwe_dev[j] += l->stat[2] < p_value_cutoff ? 1 : 0;
            }
        }
    }

    return 0;
}

void
BatchLocusProcessor::report_locus_overlap(size_t& n_sites, size_t& n_multiple_sites, ofstream *log_fh)
{
    n_sites = 0;
    n_multiple_sites = 0;

    //
    // Create a vector of maps, to record for each population, which sites overlap between different loci.
    //
    map<uint, set<uint>> final_map;
    //vector<map<uint, set<uint>>> sites (this->_mpopi->pops().size());
    for (uint pop = 0; pop < this->_mpopi->pops().size(); pop++) {
        for (const LocBin* locbin : this->_loci) {
            const CSLocus* loc = locbin->cloc;
            const LocSum* lsum = locbin->s->per_pop(pop);
            for (uint k = 0; k < loc->len; k++)
                if (lsum->nucs[k].num_indv > 0) {
                    final_map[lsum->nucs[k].bp].insert(loc->id);
                    //sites[pop][lsum->nucs[k].bp].insert(loc->id);
                }
        }
    }

    n_sites = final_map.size();

    for (auto bp : final_map) {
        if (bp.second.size() > 1) {
            n_multiple_sites++;
            if (log_fh != NULL) {
                *log_fh << "multilocus_pos"
                        << '\t' << this->_loci[0]->cloc->loc.chr() << ':' << bp.first + 1
                        << "\tloci=";
                join(bp.second, ',', *log_fh);
                *log_fh << '\n';
            }
        }
    }
}

void
LocusFilter::erase_snp(CSLocus *cloc, Datum **d, size_t n_samples, size_t snp_index)
{
    //
    // N.B. For this to work as expected we must be iterating over cloc->snps
    // in reverse, e.g. for(size_t i=cloc->snps.size(); i!=0;) {--i; ...}
    // because we're making vector::erase() calls.
    //
    // Prune out the given SNP. We need to update:
    // * `Locus::snps`
    // * `Locus::alleles`
    // * `Datum::model`
    // * `Datum::snpdata`
    // * `Datum::obshap`
    //

    //
    // Update the Datums.
    //
    uint col = cloc->snps[snp_index]->col;
    for (size_t s=0; s<n_samples; ++s) {
        if (d[s] == NULL)
            continue;

        assert(d[s]->model);
        assert(!d[s]->obshap.empty() && strlen(d[s]->obshap[0]) == cloc->snps.size());
        assert(d[s]->snpdata.size() == cloc->snps.size());
        
        //
        // Correct the model calls.
        //
        char& m = d[s]->model[col];
        assert(m == 'O' || m == 'E' || m == 'U');
        if (m == 'E' || (m == 'O' && d[s]->obshap[0][snp_index] != cloc->con[col]))
            m = 'U';
        //
        // Erase the nucleotide in the haplotypes.
        //
        for (char* hap : d[s]->obshap) {
            hap += snp_index;
            while(*hap != '\0') {
                *hap = *(hap+1);
                ++hap;
            }
        }
        //
        // Splice the depths.
        //
        d[s]->snpdata.erase(d[s]->snpdata.begin() + snp_index);
    }

    //
    // Update the CSLocus.
    //
    delete cloc->snps[snp_index];
    cloc->snps.erase(cloc->snps.begin() + snp_index);
    // Splice the haplotypes. (Why does this have to be a map?!)
    vector<pair<string,int>> alleles (cloc->alleles.begin(), cloc->alleles.end());
    cloc->alleles.clear();
    for (pair<string,int>& allele : alleles)
        allele.first.erase(snp_index, 1);
    cloc->alleles.insert(std::make_move_iterator(alleles.begin()), std::make_move_iterator(alleles.end()));
}

void
LocusFilter::batch_clear()
{
    this->_batch_total_loci    = 0;
    this->_batch_filtered_loci = 0;
    this->_batch_seen_loci     = 0;
}

void
LocusFilter::locus_seen()
{
    this->_seen_loci++;
    this->_batch_seen_loci++;
}

void
LocusFilter::locus_unsee()
{
    this->_seen_loci--;
    this->_batch_seen_loci--;
}

void
LocusFilter::keep_locus(LocBin *loc)
{
    this->_total_loci++;
    this->_batch_total_loci++;

    //
    // Count up the number of sites and variable sites.
    //
    const LocTally *t = loc->s->meta_pop();

    for (uint i = 0; i < loc->cloc->len; i++) {
        if (t->nucs[i].fixed == false)
            this->_variant_sites++;
        this->_total_sites++;
    }
}

bool
LocusFilter::whitelist_filter(size_t locus_id)
{
    if (this->_whitelist.empty() || this->_whitelist.count(locus_id)) {
        return false;
    } else {
        this->_filtered_loci++;
        this->_batch_filtered_loci++;
        return true;
    }
}

bool
LocusFilter::blacklist_filter(size_t locus_id)
{
    if (this->_blacklist.count(locus_id)) {
        this->_filtered_loci++;
        this->_batch_filtered_loci++;
        return true;
    } else {
        return false;
    }
}

void
LocusFilter::whitelist_snp_filter(LocBin& loc) const
{
    if (this->_whitelist.empty())
        return;
    CSLocus* cloc = loc.cloc;
    assert(this->_whitelist.count(cloc->id));
    const set<int>& snp_wl = this->_whitelist.at(cloc->id);
    if (snp_wl.empty())
        // Accept any SNP.
        return;
    for(size_t i=cloc->snps.size(); i!=0;) {
        --i;
        if (!snp_wl.count(cloc->snps[i]->col))
            erase_snp(cloc, loc.d, loc.sample_cnt, i);
    }
}

bool
LocusFilter::apply_filters_stacks(LocBin& loc, ostream& log_fh, const MetaPopInfo& mpopi)
{
    //
    // Apply the -r/-p thresholds at the locus level.
    //
    if (this->filter(&mpopi, loc.d, loc.cloc))
        return true;

    //
    // Remove polyallelic SNPs.
    //
    for (uint snp_index = loc.cloc->snps.size(); snp_index > 0;) {
        --snp_index;
        const SNP& snp = *loc.cloc->snps[snp_index];
        if (snp.rank_3 != 0) {
            this->_filtered_sites++;
            if (verbose) {
                log_fh << "pruned_polymorphic_site\t"
                       << loc.cloc->id << "\t"
                       << loc.cloc->loc.chr() << "\t"
                       << loc.cloc->sort_bp(snp.col) +1 << "\t"
                       << snp.col << "\t"
                       << "polyallelic\n";
            }
            this->erase_snp(loc.cloc, loc.d, mpopi.n_samples(), snp_index);
        }
    }

    //
    // Filter genotypes for depth.
    //
    this->gt_depth_filter(loc.d, loc.cloc);

    //
    // Identify individual SNPs that are below the -r threshold or the minor allele
    // frequency threshold (-a). In these cases we will remove the SNP, but keep the locus.
    //
    loc.s = new LocPopSum(strlen(loc.cloc->con), mpopi);
    this->filter_snps(loc, mpopi, log_fh);

    //
    // If write_single_snp, write_random_snp or filter_haplotype_wise has been specified,
    // mark sites to be pruned using the whitelist.
    //
    assert(write_single_snp + write_random_snp + filter_haplotype_wise <= 1);
    if (write_single_snp)
        this->keep_single_snp(loc.cloc, loc.d, mpopi.n_samples(), loc.s->meta_pop());
    else if (write_random_snp)
        this->keep_random_snp(loc.cloc, loc.d, mpopi.n_samples(), loc.s->meta_pop());
    else if (filter_haplotype_wise)
        this->filter_haps(loc, mpopi, log_fh);

    //
    // Regenerate summary statistics after pruning SNPs.
    //
    loc.s->sum_pops(loc.cloc, loc.d, mpopi, verbose, cout);
    loc.s->tally_metapop(loc.cloc);
    return false;
}

bool
LocusFilter::apply_filters_external(LocBin& loc, ostream& log_fh, const MetaPopInfo& mpopi)
{
    //
    // Apply the -r/-p thresholds at the locus level.
    //
    if (this->filter(&mpopi, loc.d, loc.cloc))
        return true;

    //
    // Identify individual SNPs that are below the -r threshold or the minor allele
    // frequency threshold (-a). In these cases we will remove the SNP, but keep the locus.
    //
    loc.s = new LocPopSum(strlen(loc.cloc->con), mpopi);
    this->filter_snps(loc, mpopi, log_fh);

    loc.s->sum_pops(loc.cloc, loc.d, mpopi, verbose, cout);
    loc.s->tally_metapop(loc.cloc);
    return false;
}

bool
LocusFilter::filter(const MetaPopInfo *mpopi, Datum **d, CSLocus *cloc)
{
    // Filter out populations that don't have enough samples.
    size_t n_pops_present = 0;
    for (const Pop& pop : mpopi->pops()) {
        size_t n_samples_present = 0;
        for (size_t s=pop.first_sample; s<=pop.last_sample; ++s)
            if (d[s] != NULL)
                ++n_samples_present;
        if ((double) n_samples_present / pop.n_samples() >= min_samples_per_pop) {
            ++n_pops_present;
        } else {
            // Remove all samples of that population.
            for (size_t s=pop.first_sample; s<=pop.last_sample; ++s) {
                if (d[s] != NULL) {
                    delete d[s];
                    d[s] = NULL;
                    cloc->cnt--;
                }
            }
        }
    }
    // Check the number of (remaining) samples.
    size_t n_samples_present = 0;
    for (size_t s=0; s<mpopi->n_samples(); ++s)
        if (d[s] != NULL)
            ++n_samples_present;
    // Determine if the locus is to be kept.
    if (n_pops_present < size_t(min_populations)
        || (double) n_samples_present / mpopi->n_samples() < min_samples_overall)
    {
        this->_filtered_loci++;
        this->_batch_filtered_loci++;
        return true;
    }
    return false;
}

void
LocusFilter::gt_depth_filter(Datum** data, const CSLocus* cloc)
{
    if (min_gt_depth == 0)
        return;
    else if (cloc->snps.empty())
        return;
    for (size_t spl=0; spl<this->_sample_cnt; ++spl) {
        Datum* d = data[spl];
        assert(d->obshap.size() == 2);
        assert(strlen(d->obshap[0]) == cloc->snps.size()
            && strlen(d->obshap[1]) == cloc->snps.size());
        for (size_t snp=0; snp<cloc->snps.size(); ++snp) {
            size_t col = cloc->snps[snp]->col;
            if (d->model[col] == 'U')
                continue;
            // Check this genotype's depth.
            if (d->snpdata[snp].tot_depth < min_gt_depth) {
                d->model[col] = 'U';
                for(size_t i=0; i<2; ++i)
                    d->obshap[i][snp] = 'N';
            }
        }
    }
}

void
LocusFilter::init(MetaPopInfo *mpopi)
{
    this->_pop_cnt    = mpopi->pops().size();
    this->_sample_cnt = mpopi->n_samples();

    assert(this->_pop_cnt > 0);
    assert(this->_sample_cnt > 0);

    if (this->_pop_order != NULL)
        delete [] this->_pop_order;
    if (this->_samples != NULL)
        delete [] this->_samples;
    if (this->_pop_tot != NULL)
        delete [] this->_pop_tot;

    this->_pop_order  = new size_t [this->_pop_cnt];
    this->_samples    = new size_t [this->_sample_cnt];
    this->_pop_tot    = new size_t [this->_pop_cnt];

    this->_filtered_loci  = 0;
    this->_total_loci     = 0;
    this->_filtered_sites = 0;
    this->_total_sites    = 0;

    size_t pop_sthg = 0;

    for (size_t i_pop = 0; i_pop < mpopi->pops().size(); ++i_pop) {
        const Pop& pop = mpopi->pops()[i_pop];
        this->_pop_tot[pop_sthg]  = 0;

        for (uint i = pop.first_sample; i <= pop.last_sample; i++) {
            this->_samples[i] = pop_sthg;
            this->_pop_tot[pop_sthg]++;
        }
        this->_pop_order[pop_sthg] = i_pop;
        pop_sthg++;
    }
}

void
LocusFilter::keep_single_snp(CSLocus* cloc, Datum** d, size_t n_samples, const LocTally* t) const
{
    //
    // Check that we have at least one variable site within this population for this locus.
    //
    size_t n_actual_snps = 0;
    for (const SNP* snp : cloc->snps)
        if (!t->nucs[snp->col].fixed)
            ++n_actual_snps;
    if (n_actual_snps == 0)
        return;

    //
    // Find the first SNP that is not fixed in this subpopulation.
    //
    size_t kept_snp = 0;
    while (kept_snp != cloc->snps.size()) {
        if (t->nucs[cloc->snps[kept_snp]->col].fixed == false)
            break;
        kept_snp++;
    }

    //
    // Remove all the SNPs except for the one marked previously.
    //
    for (size_t snp=cloc->snps.size(); snp!=0;) {
        --snp;
        if (snp != kept_snp)
            erase_snp(cloc, d, n_samples, snp);
    }
}

void
LocusFilter::keep_random_snp(CSLocus* cloc, Datum** d, size_t n_samples, const LocTally* t) const
{
    //
    // Check that we have at least one variable site within this population for this locus.
    //
    size_t n_actual_snps = 0;
    for (const SNP* snp : cloc->snps)
        if (!t->nucs[snp->col].fixed)
            ++n_actual_snps;
    if (n_actual_snps == 0)
        return;

    //
    // Identify a random SNP that isn't fixed in this subset of populations.
    //
    size_t kept_snp_i;
    do {
        kept_snp_i = rand() % cloc->snps.size();
    } while (t->nucs[cloc->snps[kept_snp_i]->col].fixed);

    //
    // Remove all the SNPs except for the one marked previously.
    //
    for (auto snp_i=cloc->snps.size(); snp_i!=0; ) {
        --snp_i;
        if (snp_i != kept_snp_i)
            erase_snp(cloc, d, n_samples, snp_i);
    }
}

void
LocusFilter::filter_snps(LocBin& loc, const MetaPopInfo& mpopi, ostream &log_fh)
{
    assert(loc.s != NULL);
    CSLocus* cloc = loc.cloc;
    Datum** d = loc.d;
    LocPopSum* s = loc.s;

    //
    // If this locus is fixed, ignore it.
    //
    if (cloc->snps.size() == 0)
        return;

    loc.s->sum_pops(loc.cloc, loc.d, mpopi, verbose, cout);
    loc.s->tally_metapop(loc.cloc);
    const LocTally *t = s->meta_pop();
    
    for (uint snp_index = cloc->snps.size(); snp_index > 0;) { // Must be reverse because we `erase()` things.
        --snp_index;
        uint col = cloc->snps[snp_index]->col;
        bool sample_prune = false;
        bool overall_sample_prune = false;
        bool maf_prune    = false;
        bool het_prune    = false;

        //
        // If the site is fixed, ignore it.
        //
        NucTally& nuct = t->nucs[col];
        if (t->nucs[col].fixed == true)
            continue;

        size_t n_pruned_pops = 0;
        size_t n_samples_left = 0;
        for (size_t p = 0; p < s->pop_cnt(); ++p) {
            const LocSum* sum = s->per_pop(p);
            if (sum->nucs[col].incompatible_site) {
                DOES_NOT_HAPPEN;
            }
            if ((double) sum->nucs[col].num_indv / this->_pop_tot[p] >= min_samples_per_pop) {
                n_samples_left += sum->nucs[col].num_indv;
            } else {
                ++n_pruned_pops;
                const Pop& pop = mpopi.pops()[p];
                for (uint k = pop.first_sample; k <= pop.last_sample; k++) {
                    if (d[k] == NULL || col >= (uint) d[k]->len)
                        continue;

                    if (d[k]->model != NULL)
                        d[k]->model[col] = 'U';
                    for(size_t j=0; j < 2; ++j)
                        d[k]->obshap[j][snp_index] = 'N';
                }
            }
        }
        if (mpopi.pops().size() - n_pruned_pops < (uint) min_populations)
            sample_prune = true;
        else if ((double) n_samples_left / mpopi.n_samples() < min_samples_overall)
            overall_sample_prune = true;

        if (t->nucs[col].allele_cnt > 1) {
            //
            // Test for minor allele frequency.
            //
            if ((1 - nuct.p_freq) < minor_allele_freq
                    || long(std::round((1 - nuct.p_freq) * 2.0 * nuct.num_indv)) < minor_allele_cnt)
                maf_prune = true;
            //
            // Test for observed heterozygosity.
            //
            if (t->nucs[col].obs_het > max_obs_het)
                het_prune = true;
        }

        if (maf_prune || het_prune || sample_prune || overall_sample_prune) {
            this->_filtered_sites++;
            if (verbose) {
                log_fh << "pruned_polymorphic_site\t"
                       << cloc->id << "\t"
                       << cloc->loc.chr() << "\t"
                       << cloc->sort_bp(col) +1 << "\t"
                       << col << "\t";
                if (sample_prune)
                    log_fh << "min_samples_per_pop\n";
                else if (overall_sample_prune)
                    log_fh << "min_samples_overall\n";
                else if (maf_prune)
                    log_fh << "maf_limit\n";
                else if (het_prune)
                    log_fh << "obshet_limit\n";
                else
                    DOES_NOT_HAPPEN;
            }
            this->erase_snp(cloc, d, mpopi.n_samples(), snp_index);
        }
    }
}

void
LocusFilter::filter_haps(LocBin& loc, const MetaPopInfo& mpopi, ostream &log_fh)
{
    CSLocus* cloc = loc.cloc;
    Datum** d = loc.d;

    // Tally SNP abundancies.
    vector<pair<size_t,size_t>> snps_n_samples;
    for (size_t i=0; i<cloc->snps.size(); ++i) {
        snps_n_samples.push_back({0, i});
        uint col = cloc->snps[i]->col;
        for (size_t s=0; s<mpopi.n_samples(); ++s)
            if (d[s] != NULL
                    && size_t(d[s]->len) > col && (d[s]->model[col] == 'O' || d[s]->model[col] == 'E'))
                ++snps_n_samples.back().first;
    }
    std::sort(snps_n_samples.rbegin(), snps_n_samples.rend());

    // Prune SNPs until enough samples have fully resolved haplotypes.
    // (N.B. This is applied after snp-wise -r/-R as any individual SNP failing
    // these filters would cause them to fail at the haplotype level as well anyway.)
    while (!cloc->snps.empty()) {
        // Tally haplotypes.
        size_t n_pops_present_hapwise = 0;
        size_t n_samples_w_haps = 0;
        for (const Pop& pop : mpopi.pops()) {
            size_t n_pop_samples = 0;
            for (size_t s=pop.first_sample; s<=pop.last_sample; ++s) {
                if (d[s] == NULL)
                    continue;
                assert(d[s]->obshap.size() == 2);
                if (strchr(d[s]->obshap[0], 'N') == NULL)
                    ++n_pop_samples;
            }
            if ((double) n_pop_samples / pop.n_samples() >= min_samples_per_pop) {
                ++n_pops_present_hapwise;
                n_samples_w_haps += n_pop_samples;
            }
        }
        // Check the number of haplotypes.
        if (n_pops_present_hapwise >= size_t(min_populations)
                && (double) n_samples_w_haps / mpopi.n_samples() >= min_samples_overall)
            // Filters are satisfied.
            break;
        // Prune the SNP that is present in the fewest samples.
        ++this->_filtered_sites;
        assert(!snps_n_samples.empty());
        size_t snp_i = snps_n_samples.back().second;
        this->erase_snp(cloc, d, mpopi.n_samples(), snp_i);
        snps_n_samples.pop_back();
        assert(snps_n_samples.size() == cloc->snps.size());
        for (pair<size_t,size_t>& snp : snps_n_samples)
            if (snp.second > snp_i)
                --snp.second;
    }
}

int
CatalogDists::accumulate_pre_filtering(const size_t sample_cnt, const CSLocus *loc)
{
    size_t missing;

    if (this->_pre_valid.count(loc->cnt) == 0)
        this->_pre_valid[loc->cnt] = 1;
    else
        this->_pre_valid[loc->cnt]++;

    missing = sample_cnt - loc->cnt;

    if (this->_pre_absent.count(missing) == 0)
        this->_pre_absent[missing] = 1;
    else
        this->_pre_absent[missing]++;

    if (this->_pre_snps_per_loc.count(loc->snps.size()) == 0)
        this->_pre_snps_per_loc[loc->snps.size()] = 1;
    else
        this->_pre_snps_per_loc[loc->snps.size()]++;

    return 0;
}

int
CatalogDists::accumulate(const vector<LocBin *> &loci)
{
    const CSLocus *loc;
    const LocTally *t;
    size_t missing;

    for (uint i = 0; i < loci.size(); i++) {
        loc = loci[i]->cloc;

        if (this->_post_valid.count(loc->cnt) == 0)
            this->_post_valid[loc->cnt] = 1;
        else
            this->_post_valid[loc->cnt]++;

        missing = loci[i]->sample_cnt - loc->cnt;

        if (this->_post_absent.count(missing) == 0)
            this->_post_absent[missing] = 1;
        else
            this->_post_absent[missing]++;

        //
        // Don't count SNPs that are fixed in the metapopulation.
        //
        t = loci[i]->s->meta_pop();
        size_t n_actual_snps = 0;
        for (const SNP* snp : loc->snps)
            if (!t->nucs[snp->col].fixed)
                ++n_actual_snps;

        if (this->_post_snps_per_loc.count(n_actual_snps) == 0)
            this->_post_snps_per_loc[n_actual_snps] = 1;
        else
            this->_post_snps_per_loc[n_actual_snps]++;
    }

    return 0;
}

int
CatalogDists::write_results(ostream &log_fh)
{
    map<size_t, size_t>::iterator cnt_it;
    string section;
    auto begin_section = [&](const string& s){
        section = s;
        log_fh << "\n" << "BEGIN " << section << "\n";
    };
    auto end_section = [&](){
        log_fh << "END " << section << "\n";
    };

    begin_section("samples_per_loc_prefilters");
    log_fh << "# Distribution of valid samples matched to a catalog locus prior to filtering.\n"
           << "n_samples\tn_loci\n";
    for (cnt_it = this->_pre_valid.begin(); cnt_it != this->_pre_valid.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";
    end_section();

    begin_section("missing_samples_per_loc_prefilters");
    log_fh << "# Distribution of missing samples for each catalog locus prior to filtering.\n"
           << "# Absent samples at locus\tCount\n";
    for (cnt_it = this->_pre_absent.begin(); cnt_it != this->_pre_absent.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";
    end_section();

    begin_section("snps_per_loc_prefilters");
    log_fh << "# Distribution of the number of SNPs per catalog locus prior to filtering.\n"
           << "n_snps\tn_loci\n";
    for (cnt_it = this->_pre_snps_per_loc.begin(); cnt_it != this->_pre_snps_per_loc.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";
    end_section();

    begin_section("samples_per_loc_postfilters");
    log_fh << "# Distribution of valid samples matched to a catalog locus after filtering.\n"
           << "n_samples\tn_loci\n";
    for (cnt_it = this->_post_valid.begin(); cnt_it != this->_post_valid.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";
    end_section();

    begin_section("missing_samples_per_loc_postfilters");
    log_fh << "# Distribution of missing samples for each catalog locus after filtering.\n"
           << "# Absent samples at locus\tCount\n";
    for (cnt_it = this->_post_absent.begin(); cnt_it != this->_post_absent.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";
    end_section();

    begin_section("snps_per_loc_postfilters");
    log_fh << "# Distribution of the number of SNPs per catalog locus (after filtering).\n"
           << "n_snps\tn_loci\n";
    for (cnt_it = this->_post_snps_per_loc.begin(); cnt_it != this->_post_snps_per_loc.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";
    end_section();

    return 0;
}

/*
int
merge_shared_cutsite_loci(map<int, CSLocus *> &catalog,
                          PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum,
                          map<int, pair<merget, int> > &merge_map,
                          ofstream &log_fh)
{
    map<string, vector<CSLocus *> >::iterator it;
    CSLocus *cur, *next;
    Datum  **d_1, **d_2;
    double   prune_pct;
    uint unmergable, tot_loci, tot_samp;
    uint success           = 0;
    uint failure           = 0;
    uint overlap           = 0;
    uint simple_merge_cnt  = 0;
    uint complex_merge_cnt = 0;
    uint missing_samps_cnt = 0;
    uint phase_fail_cnt    = 0;
    uint nomapping_cnt     = 0;
    uint multimapping_cnt  = 0;
    uint multifails_cnt    = 0;

    tot_loci = pmap->loci_cnt();

    set<int> loci_to_destroy;
    map<int, int> missing_samps_dist;

    cout << "To merge adjacent loci at least " << merge_prune_lim * 100 << "% of samples must have both adjacent loci;"
         << " the remaining " << 100 - (merge_prune_lim * 100) << "% of individuals will be pruned.\n"
         << "Attempting to merge adjacent loci that share a cutsite...";

    if (verbose)
        log_fh << "\n#\n# List of locus pairs that share a cutsite that failed to merge because they could not be phased.\n#\n";

    //
    // Iterate over each chromosome.
    //
    for (it = pmap->ordered_loci_nconst().begin(); it != pmap->ordered_loci_nconst().end(); it++) {
        //
        // Iterate over each ordered locus on this chromosome.
        //
        next = it->second[0];
        for (uint pos = 1; pos < it->second.size(); pos++) {
            cur  = next;
            next = it->second[pos];

            //
            // Do these two loci overlap?
            //   +Must occur on opposite strands
            //   +Must overlap according to the length of the cutsite.
            //
            if (((cur->loc.strand == strand_minus && next->loc.strand == strand_plus) &&
                 ((int) (cur->loc.bp  - next->loc.bp + 1) == renz_olap[enz])) ||
                ((cur->loc.strand == strand_plus  && next->loc.strand == strand_minus) &&
                 ((int) (next->loc.bp - cur->loc.bp  + 1) == renz_olap[enz]))) {
                overlap++;

                d_1        = pmap->locus(cur->id);
                d_2        = pmap->locus(next->id);
                unmergable = 0;
                tot_samp   = 0;

                //
                // Check if all members of the population contain these two loci (or are missing both).
                //
                for (int i = 0; i < pmap->sample_cnt(); i++) {
                    if (d_1[i] != NULL || d_2[i] != NULL)
                        tot_samp++;
                    if ((d_1[i] != NULL && d_2[i] == NULL) ||
                        (d_1[i] == NULL && d_2[i] != NULL))
                        unmergable++;
                }

                prune_pct = (double) (tot_samp - unmergable) / (double) tot_samp;

                //
                // If some of the individuals only have one locus and not the other, prune them out.
                //
                if (prune_pct < 1.0 && prune_pct >= merge_prune_lim) {
                    for (int i = 0; i < pmap->sample_cnt(); i++)
                        if (d_1[i] != NULL && d_2[i] == NULL) {
                            delete d_1[i];
                            d_1[i] = NULL;
                        } else if (d_1[i] == NULL && d_2[i] != NULL) {
                            delete d_2[i];
                            d_2[i] = NULL;
                        }
                }

                //
                // If possible, merge the two loci together.
                //
                if (prune_pct < merge_prune_lim) {
                    int pct = (int) (prune_pct * 100);
                    missing_samps_dist[pct]++;
                    if (verbose) log_fh << "Missing samples, Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; "
                                        << pct << "% present (" << 100 - pct << "% missing)\n";
                    missing_samps_cnt++;
                    failure++;
                    continue;
                }

                phaset res = merge_and_phase_loci(pmap, cur, next, loci_to_destroy, log_fh);
                switch(res) {
                case multiple_fails:
                    if (verbose) log_fh << "Failed to phase, Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; "
                                        << "multiple failures\n";
                    multifails_cnt++;
                    phase_fail_cnt++;
                    failure++;
                    break;
                case multimapping_fail:
                    if (verbose) log_fh << "Failed to phase, Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; "
                                        << "multimapping in one or more individuals\n";
                    multimapping_cnt++;
                    phase_fail_cnt++;
                    failure++;
                    break;
                case nomapping_fail:
                    if (verbose) log_fh << "Failed to phase, Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; "
                                        << "no mapping in one or more individuals\n";
                    nomapping_cnt++;
                    phase_fail_cnt++;
                    failure++;
                    break;
                case complex_phase:
                    if (verbose) log_fh << "Phased Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; "
                                        << "a complex phasing operation.\n";
                    complex_merge_cnt++;
                    success++;
                    merge_map[cur->id]  = make_pair(merge_sink, next->id);
                    merge_map[next->id] = make_pair(merge_src, cur->id);
                    break;
                case simple_merge:
                    if (verbose) log_fh << "Phased Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; "
                                        << "a simple merge operation.\n";
                    simple_merge_cnt++;
                    success++;
                    merge_map[cur->id]  = make_pair(merge_sink, next->id);
                    merge_map[next->id] = make_pair(merge_src, cur->id);
                    break;
                default:
                    cerr << "Warning: Merge failure.\n";
                    break;
                }
            }
        }
    }

    //
    // Remove those loci that have been merged from both the popualtion map and catalog.
    //
    set<int> emptyset;
    pmap->prune(loci_to_destroy);
    reduce_catalog(catalog, emptyset, loci_to_destroy);

    cout << "done.\n"
         << "Of " << tot_loci << " loci, "
         << overlap << " pairs share a cutsite; "
         << success << " pairs were merged; "
         << failure << " pairs failed to merge; "
         << pmap->loci_cnt() << " remaining loci.\n"
         << "  Of those merged, " << simple_merge_cnt << " required only a simple merge without phasing; "
         << "while " << complex_merge_cnt << " required phasing.\n"
         << "  Of those that failed to merge, " << missing_samps_cnt << " were missing one of the two haplotypes in one or more samples; "
         << "while " << phase_fail_cnt << " failed to be phased.\n"
         << "    Of those that failed to phase, " << nomapping_cnt << " failed due to a lack of haplotype mappings; "
         << multimapping_cnt << " failed due to multiple haplotype mappings; " << multifails_cnt << " failed due to both.\n";

    log_fh << "\n#\n# Merging adjacent loci with a shared restriction enzyme cutsite\n#\n"
           << "Of " << tot_loci << " loci, "
           << overlap << " pairs share a cutsite; "
           << success << " pairs were merged; "
           << failure << " pairs failed to merge; "
           << pmap->loci_cnt() << " remaining loci.\n"
           << "  Of those merged, " << simple_merge_cnt << " required only a simple merge without phasing; "
           << "while " << complex_merge_cnt << " required phasing.\n"
           << "  Of those that failed to merge, " << missing_samps_cnt << " were missing one of the two haplotypes in one or more samples; "
           << "while " << phase_fail_cnt << " failed to be phased.\n"
           << "    Of those that failed to phase, " << nomapping_cnt << " failed due to a lack of haplotype mappings; "
           << multimapping_cnt << " failed due to multiple haplotype mappings; " << multifails_cnt << " failed due to both.\n";
    log_fh << "#\n# Distribution of loci with samples missing one of two loci to be merged\n"
           << "# Percent samples with both loci present\tNumber of cases\n";
    map<int, int>::iterator mit;
    for (mit = missing_samps_dist.begin(); mit != missing_samps_dist.end(); mit++)
        log_fh << mit->first << "\t" << mit->second << "\n";
    log_fh << "\n";

    return 0;
}
*/
/*phaset
merge_and_phase_loci(PopMap<CSLocus> *pmap, CSLocus *cur, CSLocus *next,
                     set<int> &loci_to_destroy,
                     ofstream &log_fh)
{
    Datum **d_1 = pmap->locus(cur->id);
    Datum **d_2 = pmap->locus(next->id);

    set<int>    phased_results;
    set<string> phased_haplotypes;
    string      merged_hap;
    char       *h_1, *h_2;
    int         merge_type;

    if (verbose) log_fh << "Attempting to phase source locus " << cur->id << " with sink locus " << next->id << "\n";

    int sample_cnt        = 0;
    int phased_sample_cnt = 0;
    //
    // Take a census of the already phased haplotypes. We have phased haplotypes
    // if for individual i:
    //   1. d_1 has a single haplotype and d_2 has a single haplotype
    //   2. d_1 has a single haplotpye and d_2 has multiple haplotypes
    //   3. d_1 has multiple haplotpyes and d_2 has a single haplotype
    //
    // If one or both of the loci have no SNPs, then the haplotype is
    // recorded as "consensus." Check that condition before we start merging.
    //
    if (cur->snps.size() > 0 && next->snps.size() > 0)
        merge_type = 0;
    else if (cur->snps.size() == 0)
        merge_type = 1;
    else if (next->snps.size() == 0)
        merge_type = 2;
    else
        merge_type = 3;

    for (int i = 0; i < pmap->sample_cnt(); i++) {
        if (d_1[i] == NULL || d_2[i] == NULL)
            continue;
        else if (d_1[i]->obshap.size() > 1 && d_2[i]->obshap.size() > 1)
            continue;
        else {
            for (uint j = 0; j < d_1[i]->obshap.size(); j++) {
                for (uint k = 0; k < d_2[i]->obshap.size(); k++) {
                    switch (merge_type) {
                    case 0:
                        merged_hap = string(d_1[i]->obshap[j]) + string(d_2[i]->obshap[k]);
                        break;
                    case 1:
                        merged_hap = string(d_2[i]->obshap[k]);
                        break;
                    case 2:
                        merged_hap = string(d_1[i]->obshap[j]);
                        break;
                    case 3:
                    default:
                        merged_hap = "consensus";
                        break;
                    }
                    phased_haplotypes.insert(merged_hap);
                    // cout << "Phasing: '" << d_1[i]->obshap[j] << "' + '" << d_2[i]->obshap[k] << "' => '" << merged_hap << "'\n";
                }
            }
            phased_sample_cnt++;
            sample_cnt++;
        }
    }

    //
    // Indicate that these two loci had a simple merge, with no phasing necessary.
    //
    phased_results.insert(simple_merge);

    //
    // Now we need to check if we can phase the remaining haplotypes.
    //
    for (int i = 0; i < pmap->sample_cnt(); i++) {
        if (d_1[i] == NULL || d_2[i] == NULL)
            continue;
        else if (d_1[i]->obshap.size() > 1 && d_2[i]->obshap.size() > 1) {
            // cout << "Attempting to phase individual " << i << ": " << d_1[i]->id << " / " << d_2[i]->id << "\n";

            sample_cnt++;
            //
            // We should be able to find a sinlge phasing mapping for each haplotype from d_1 to d_2
            // that includes all the haplotypes in these two loci.
            //
            vector<pair<char *, char *> > seen_phased;
            uint tot_obshap = d_1[i]->obshap.size() + d_2[i]->obshap.size();
            uint phased_cnt = 0;
            for (uint j = 0; j < d_1[i]->obshap.size(); j++) {
                for (uint k = 0; k < d_2[i]->obshap.size(); k++) {
                    // cout << "  " << d_1[i]->obshap[j] << " + " << d_2[i]->obshap[k];
                    //
                    // Record each pair of haplotypes that has been seen phased previously.
                    //
                    if (phased_haplotypes.count(string(d_1[i]->obshap[j]) + string(d_2[i]->obshap[k]))) {
                        seen_phased.push_back(make_pair(d_1[i]->obshap[j], d_2[i]->obshap[k]));
                        // cout << " => " << d_1[i]->obshap[j] << d_2[i]->obshap[k];
                    }
                    // cout << "\n";
                }
            }
            //
            // Now, we will iterate over all sets of phased haplotypes and look
            // for combinations that use all four individual haplotypes.
            //
            for (uint j = 0; j < seen_phased.size(); j++) {
                for (uint k = j; k < seen_phased.size(); k++) {
                    set<char *> incorporated_haplotypes;
                    //
                    // Count the number of distinct char pointers. If this combination
                    // of haplotypes includes all unphased haplotypes, count it.
                    //
                    incorporated_haplotypes.insert(seen_phased[j].first);
                    incorporated_haplotypes.insert(seen_phased[j].second);
                    incorporated_haplotypes.insert(seen_phased[k].first);
                    incorporated_haplotypes.insert(seen_phased[k].second);
                    if (incorporated_haplotypes.size() == tot_obshap)
                        phased_cnt++;
                }
            }

            //
            // If one pair of haplotypes is mapped, but the other is not, assume the second pair or
            // haplotypes must be phased by process of elimination.
            //
            if (phased_cnt == 0 && seen_phased.size() == 1) {
                h_1 = seen_phased[0].first  == d_1[i]->obshap[1] ?
                    d_1[i]->obshap[0] : d_1[i]->obshap[1];
                h_2 = seen_phased[0].second == d_2[i]->obshap[1] ?
                    d_2[i]->obshap[0] : d_2[i]->obshap[1];
                phased_haplotypes.insert(string(h_1) + string(h_2));
                phased_cnt++;
                // cout << "  Phasing: '" << hap_1 << "' + '" << hap_2 << "' => '" << string(hap_1) + string(hap_2) << "'\n";
            }

            if (phased_cnt == 0) {
                phased_results.insert(nomapping_fail);
                if (verbose) log_fh << "    Locus NOT phased in individual " << i << "; loci " << d_1[i]->id << " / " << d_2[i]->id << " no mapping found.\n";
            } else if (phased_cnt == 1) {
                phased_sample_cnt++;
                phased_results.insert(complex_phase);
            } else {
                phased_results.insert(multimapping_fail);
                if (verbose) log_fh << "    Locus NOT phased in individual " << i << "; loci " << d_1[i]->id << " / " << d_2[i]->id << " multiple mappings found.\n";
            }
        }
    }

    if (phased_sample_cnt != sample_cnt) {
        if (phased_results.count(nomapping_fail) > 0 &&
            phased_results.count(multimapping_fail) > 0)
            return multiple_fails;
        else if (phased_results.count(nomapping_fail) > 0)
            return nomapping_fail;
        else if (phased_results.count(multimapping_fail) > 0)
            return multimapping_fail;
        else {
            cout << "WE SHOULD NOT GET HERE\n";
            return merge_failure;
        }
    }

    //
    // Okay, merge these two loci together.
    //
    if (!merge_datums(pmap->sample_cnt(), cur->len, d_1, d_2, phased_haplotypes, merge_type))
        return merge_failure;

    //
    // Merge the catalog entries together.
    //
    if (!merge_csloci(cur, next, phased_haplotypes))
        return merge_failure;

    //
    // Mark the merged locus for destruction.
    //
    loci_to_destroy.insert(next->id);

    if (phased_results.count(complex_phase) > 0)
        return complex_phase;
    return simple_merge;
}*/

/*int
merge_csloci(CSLocus *sink, CSLocus *src, set<string> &phased_haplotypes)
{
    //
    // We assume that we are merging two loci: one on the negative strand, one on the
    // positive. We will keep the sink cslocus and delete the src cslocus.
    //   -> The sink cslocus is assumed to be on the negative strand.
    //

    //
    // 1. Reverse complement the SNP coordinates in the sink locus so that they are
    //    enumerated on the positive strand. Complement the alleles as well.
    //
    for (uint j = 0; j < sink->snps.size(); j++) {
        sink->snps[j]->col    = sink->len - sink->snps[j]->col - 1;
        sink->snps[j]->rank_1 = reverse(sink->snps[j]->rank_1);
        sink->snps[j]->rank_2 = reverse(sink->snps[j]->rank_2);
        sink->snps[j]->rank_3 = reverse(sink->snps[j]->rank_3);
        sink->snps[j]->rank_4 = reverse(sink->snps[j]->rank_4);
    }

    //
    // 2. Adjust the SNP coordinates in the src locus to account for the now, longer length.
    //
    for (uint j = 0; j < src->snps.size(); j++)
        src->snps[j]->col = sink->len + src->snps[j]->col - renz_olap[enz];

    //
    // 3. Combine SNPs between the two catalog loci: add the SNPs from the sink (formerly on the
    //    negative strand) in reverse order, followed by the SNPs from the src.
    //
    vector<SNP *> tmpsnp;
    for (int j = (int) sink->snps.size() - 1; j >= 0; j--)
        tmpsnp.push_back(sink->snps[j]);
    for (uint j = 0; j < src->snps.size(); j++)
        tmpsnp.push_back(src->snps[j]);
    sink->snps.clear();
    for (uint j = 0; j < tmpsnp.size(); j++)
        sink->snps.push_back(tmpsnp[j]);

    //
    // 4. Adjust the genomic location of the sink locus.
    //
    uint bp = sink->sort_bp();
    sink->loc.bp     = bp;
    sink->loc.strand = strand_plus;

    //
    // 5. Adjust the length of the sequence.
    //
    sink->len += src->len - renz_olap[enz];

    //
    // 6. Merge the consensus sequence together.
    //
    char *new_con = rev_comp(sink->con);
    delete [] sink->con;
    sink->con = new_con;
    new_con   = new char[sink->len + 1];
    strcpy(new_con, sink->con);
    delete [] sink->con;
    sink->con = new_con;
    new_con  += src->len - renz_olap[enz];
    strcpy(new_con, src->con);

    //
    // 7. Record the now phased haplotypes.
    //
    sink->alleles.clear();
    set<string>::iterator it;
    for (it = phased_haplotypes.begin(); it != phased_haplotypes.end(); it++)
        sink->alleles[*it] = 0;

    // cout << "CSLocus " << sink->id << ":\n"
    //   << "Length: " << sink->len << "; Chr: " << sink->loc.chr << "; BP: " << sink->sort_bp() << "; strand: " << (sink->loc.strand == strand_plus ? "+" : "-") << "\n"
    //   << "  SNPs:\n";
    // for (uint j = 0; j < sink->snps.size(); j++)
    //  cout << "    Col: " << sink->snps[j]->col
    //       << "    Rank 1: " << sink->snps[j]->rank_1
    //       << "    Rank 2: " << sink->snps[j]->rank_2 << "\n";
    // cout << "  Alleles:\n";
    // map<string, int>::iterator ait;
    // for (ait = sink->alleles.begin(); ait != sink->alleles.end(); ait++)
    //  cout << "    " << ait->first << "\n";

    return 1;
}*/

/*int
merge_datums(int sample_cnt,
             int sink_locus_len,
             Datum **sink, Datum **src,
             set<string> &phased_haplotypes,
             int merge_type)
{
    char           tmphap[id_len], *new_hap;
    uint           haplen, model_len, offset;
    vector<SNP *>  tmpsnp;
    vector<string> tmpobshap;
    vector<int>    tmpobsdep;

    //
    // We assume that we are merging two loci: one on the negative strand, one on the
    // positive. We will keep the sink datum and delete the src datum.
    //   -The sink datum is assumed to be on the negative strand.
    //
    for (int i = 0; i < sample_cnt; i++) {
        if (sink[i] == NULL && src[i] == NULL)
            continue;
        else if (sink[i] == NULL || src[i] == NULL)
            cout << "Unexpected condition in merging datums: one datum is NULL while the other is not.\n";

        //
        // 1. Reverse complement the observed haplotypes in the sink locus.
        //
        haplen = strlen(sink[i]->obshap[0]);
        for (uint j = 0; j < sink[i]->obshap.size(); j++) {
            for (uint k = 0; k < haplen; k++)
                tmphap[k] = reverse(sink[i]->obshap[j][haplen - k - 1]);
            tmphap[haplen] = '\0';
            strcpy(sink[i]->obshap[j], tmphap);
        }
    }

    //
    // 2. Combine observed haplotypes between the two datums while phasing them.
    //    2.1 First combine the haplotypes from samples that are already in phase.
    //
    string      merged_hap;
    vector<int> to_be_phased;
    phased_haplotypes.clear();
    for (int i = 0; i < sample_cnt; i++) {
        if (sink[i] == NULL && src[i] == NULL)
            continue;

        if (sink[i]->obshap.size() > 1 && src[i]->obshap.size() > 1) {
            to_be_phased.push_back(i);
            continue;
        } else {
            tmpobshap.clear();
            tmpobsdep.clear();
            for (uint j = 0; j < sink[i]->obshap.size(); j++) {
                for (uint k = 0; k < src[i]->obshap.size(); k++) {
                    switch (merge_type) {
                    case 0:
                        merged_hap = string(sink[i]->obshap[j]) + string(src[i]->obshap[k]);
                        break;
                    case 1:
                        merged_hap = string(src[i]->obshap[j]);
                        break;
                    case 2:
                        merged_hap = string(sink[i]->obshap[j]);
                        break;
                    case 3:
                    default:
                        merged_hap = "consensus";
                        break;
                    }
                    phased_haplotypes.insert(merged_hap);
                    tmpobshap.push_back(merged_hap);
                    tmpobsdep.push_back((sink[i]->depth[j] + src[i]->depth[k]) / 2);
                }
            }
            sink[i]->depth.clear();
            for (uint j = 0; j < sink[i]->obshap.size(); j++)
                delete [] sink[i]->obshap[j];
            sink[i]->obshap.clear();
            for (uint j = 0; j < tmpobshap.size(); j++) {
                new_hap = new char[tmpobshap[j].length() + 1];
                strcpy(new_hap, tmpobshap[j].c_str());
                sink[i]->obshap.push_back(new_hap);
                sink[i]->depth.push_back(tmpobsdep[j]);
            }
        }
    }
    //
    //    2.2 Phase and combine the haplotypes from the remaining samples.
    //
    int index;
    for (uint i = 0; i < to_be_phased.size(); i++) {
        index = to_be_phased[i];
        tmpobshap.clear();
        tmpobsdep.clear();

        vector<pair<char *, char *> > seen_phased;
        uint tot_obshap = sink[index]->obshap.size() + src[index]->obshap.size();

        for (uint j = 0; j < sink[index]->obshap.size(); j++) {
            for (uint k = 0; k < src[index]->obshap.size(); k++) {
                if (phased_haplotypes.count(string(sink[index]->obshap[j]) + string(src[index]->obshap[k])))
                    seen_phased.push_back(make_pair(sink[index]->obshap[j], src[index]->obshap[k]));
            }
        }

        for (uint j = 0; j < seen_phased.size(); j++) {
            for (uint k = j; k < seen_phased.size(); k++) {
                set<char *> incorporated_haplotypes;
                incorporated_haplotypes.insert(seen_phased[j].first);
                incorporated_haplotypes.insert(seen_phased[j].second);
                incorporated_haplotypes.insert(seen_phased[k].first);
                incorporated_haplotypes.insert(seen_phased[k].second);
                if (incorporated_haplotypes.size() == tot_obshap) {
                    tmpobshap.push_back(string(seen_phased[j].first) + string(seen_phased[j].second));
                    tmpobshap.push_back(string(seen_phased[k].first) + string(seen_phased[k].second));
                    //tmpobsdep.push_back((sink[index]->depth[j] + src[index]->depth[k]) / 2);
                }
            }
        }

        sink[index]->depth.clear();
        for (uint j = 0; j < sink[index]->obshap.size(); j++)
            delete [] sink[index]->obshap[j];
        sink[index]->obshap.clear();
        for (uint j = 0; j < tmpobshap.size(); j++) {
            new_hap = new char[tmpobshap[j].length() + 1];
            strcpy(new_hap, tmpobshap[j].c_str());
            sink[index]->obshap.push_back(new_hap);
            // sink[index]->depth.push_back(tmpobsdep[j]);
        }
    }

    //
    // 3. Merge model calls; Set the length; combine the two depth and lnl measures together.
    //
    string model_calls;
    char  *p;

    for (int i = 0; i < sample_cnt; i++) {
        if (sink[i] == NULL && src[i] == NULL)
            continue;

        //
        // Merge the two strings of model calls together.
        // We need to check if the locus for this individual is shorter than the catalog
        // locus. If so, we need to expand out the model call array to be the proper length.
        //
        reverse_string(sink[i]->model);
        offset = 0;
        model_calls.clear();
        if (sink_locus_len > sink[i]->len) {
            offset = sink_locus_len - sink[i]->len;
            model_calls.assign(offset, 'N');
        }
        model_len = offset + sink[i]->len + src[i]->len - renz_olap[enz];
        model_calls.append(sink[i]->model);
        delete [] sink[i]->model;
        sink[i]->model = new char[model_len + 1];
        strcpy(sink[i]->model, model_calls.c_str());
        p  = sink[i]->model;
        p += offset + sink[i]->len - renz_olap[enz];
        strcpy(p, src[i]->model);

        sink[i]->len       = model_len;
        sink[i]->tot_depth = (sink[i]->tot_depth + src[i]->tot_depth) / 2;
        sink[i]->lnl       = (sink[i]->lnl + src[i]->lnl) / 2.0;

        //
        // Record which datum was merged into this one.
        //
        sink[i]->merge_partner = src[i]->id;
    }

    return 1;
}*/

SumStatsSummary::SumStatsSummary(size_t pop_cnt)
{
    this->_pop_cnt           = pop_cnt;
    this->_private_cnt       = new size_t[this->_pop_cnt];
    this->_sig_hwe_dev       = new size_t[this->_pop_cnt];
    this->_n                 = new double[this->_pop_cnt];
    this->_var_sites         = new double[this->_pop_cnt];

    this->_num_indv_mean     = new double[this->_pop_cnt];
    this->_num_indv_acc_mean = new double[this->_pop_cnt];
    this->_num_indv_var      = new double[this->_pop_cnt];
    this->_p_mean            = new double[this->_pop_cnt];
    this->_p_acc_mean        = new double[this->_pop_cnt];
    this->_p_var             = new double[this->_pop_cnt];
    this->_obs_het_mean      = new double[this->_pop_cnt];
    this->_obs_het_acc_mean  = new double[this->_pop_cnt];
    this->_obs_het_var       = new double[this->_pop_cnt];
    this->_obs_hom_mean      = new double[this->_pop_cnt];
    this->_obs_hom_acc_mean  = new double[this->_pop_cnt];
    this->_obs_hom_var       = new double[this->_pop_cnt];
    this->_exp_het_mean      = new double[this->_pop_cnt];
    this->_exp_het_acc_mean  = new double[this->_pop_cnt];
    this->_exp_het_var       = new double[this->_pop_cnt];
    this->_exp_hom_mean      = new double[this->_pop_cnt];
    this->_exp_hom_acc_mean  = new double[this->_pop_cnt];
    this->_exp_hom_var       = new double[this->_pop_cnt];
    this->_pi_mean           = new double[this->_pop_cnt];
    this->_pi_acc_mean       = new double[this->_pop_cnt];
    this->_pi_var            = new double[this->_pop_cnt];
    this->_fis_mean          = new double[this->_pop_cnt];
    this->_fis_acc_mean      = new double[this->_pop_cnt];
    this->_fis_var           = new double[this->_pop_cnt];

    this->_n_all                 = new double[this->_pop_cnt];
    this->_num_indv_mean_all     = new double[this->_pop_cnt];
    this->_num_indv_acc_mean_all = new double[this->_pop_cnt];
    this->_num_indv_var_all      = new double[this->_pop_cnt];
    this->_p_mean_all            = new double[this->_pop_cnt];
    this->_p_acc_mean_all        = new double[this->_pop_cnt];
    this->_p_var_all             = new double[this->_pop_cnt];
    this->_obs_het_mean_all      = new double[this->_pop_cnt];
    this->_obs_het_acc_mean_all  = new double[this->_pop_cnt];
    this->_obs_het_var_all       = new double[this->_pop_cnt];
    this->_obs_hom_mean_all      = new double[this->_pop_cnt];
    this->_obs_hom_acc_mean_all  = new double[this->_pop_cnt];
    this->_obs_hom_var_all       = new double[this->_pop_cnt];
    this->_exp_het_mean_all      = new double[this->_pop_cnt];
    this->_exp_het_acc_mean_all  = new double[this->_pop_cnt];
    this->_exp_het_var_all       = new double[this->_pop_cnt];
    this->_exp_hom_mean_all      = new double[this->_pop_cnt];
    this->_exp_hom_acc_mean_all  = new double[this->_pop_cnt];
    this->_exp_hom_var_all       = new double[this->_pop_cnt];
    this->_pi_mean_all           = new double[this->_pop_cnt];
    this->_pi_acc_mean_all       = new double[this->_pop_cnt];
    this->_pi_var_all            = new double[this->_pop_cnt];
    this->_fis_mean_all          = new double[this->_pop_cnt];
    this->_fis_acc_mean_all      = new double[this->_pop_cnt];
    this->_fis_var_all           = new double[this->_pop_cnt];

    this->_sq_n     = new double[this->_pop_cnt];
    this->_sq_n_all = new double[this->_pop_cnt];

    for (uint j = 0; j < this->_pop_cnt; j++) {
        this->_private_cnt[j]       = 0;
        this->_sig_hwe_dev[j]       = 0;
        this->_n[j]                 = 0.0;
        this->_var_sites[j]         = 0.0;
        this->_num_indv_mean[j]     = 0.0;
        this->_num_indv_acc_mean[j] = 0.0;
        this->_num_indv_var[j]      = 0.0;
        this->_p_mean[j]            = 0.0;
        this->_p_acc_mean[j]        = 0.0;
        this->_p_var[j]             = 0.0;
        this->_obs_het_mean[j]      = 0.0;
        this->_obs_het_acc_mean[j]  = 0.0;
        this->_obs_het_var[j]       = 0.0;
        this->_obs_hom_mean[j]      = 0.0;
        this->_obs_hom_acc_mean[j]  = 0.0;
        this->_obs_hom_var[j]       = 0.0;
        this->_exp_het_mean[j]      = 0.0;
        this->_exp_het_acc_mean[j]  = 0.0;
        this->_exp_het_var[j]       = 0.0;
        this->_exp_hom_mean[j]      = 0.0;
        this->_exp_hom_acc_mean[j]  = 0.0;
        this->_exp_hom_var[j]       = 0.0;
        this->_pi_mean[j]           = 0.0;
        this->_pi_acc_mean[j]       = 0.0;
        this->_pi_var[j]            = 0.0;
        this->_fis_mean[j]          = 0.0;
        this->_fis_acc_mean[j]      = 0.0;
        this->_fis_var[j]           = 0.0;

        this->_n_all[j]                 = 0.0;
        this->_num_indv_mean_all[j]     = 0.0;
        this->_num_indv_acc_mean_all[j] = 0.0;
        this->_num_indv_var_all[j]      = 0.0;
        this->_p_mean_all[j]            = 0.0;
        this->_p_acc_mean_all[j]        = 0.0;
        this->_p_var_all[j]             = 0.0;
        this->_obs_het_mean_all[j]      = 0.0;
        this->_obs_het_acc_mean_all[j]  = 0.0;
        this->_obs_het_var_all[j]       = 0.0;
        this->_obs_hom_mean_all[j]      = 0.0;
        this->_obs_hom_acc_mean_all[j]  = 0.0;
        this->_obs_hom_var_all[j]       = 0.0;
        this->_exp_het_mean_all[j]      = 0.0;
        this->_exp_het_acc_mean_all[j]  = 0.0;
        this->_exp_het_var_all[j]       = 0.0;
        this->_exp_hom_mean_all[j]      = 0.0;
        this->_exp_hom_acc_mean_all[j]  = 0.0;
        this->_exp_hom_var_all[j]       = 0.0;
        this->_pi_mean_all[j]           = 0.0;
        this->_pi_acc_mean_all[j]       = 0.0;
        this->_pi_var_all[j]            = 0.0;
        this->_fis_mean_all[j]          = 0.0;
        this->_fis_acc_mean_all[j]      = 0.0;
        this->_fis_var_all[j]           = 0.0;

        this->_sq_n[j]     = 0.0;
        this->_sq_n_all[j] = 0.0;
    }

    this->_locus_n         = 0.0;
    this->_locus_overlap_n = 0.0;
    this->_locus_pe_ctg_n  = 0.0;
    this->_locus_len_mean     = 0.0;
    this->_locus_len_acc_mean = 0.0;
    this->_locus_len_var      = 0.0;
    this->_overlap_mean     = 0.0;
    this->_overlap_acc_mean = 0.0;
    this->_overlap_var      = 0.0;
    this->_locus_gt_sites_mean     = 0.0;
    this->_locus_gt_sites_acc_mean = 0.0;
    this->_locus_gt_sites_var      = 0.0;
}

int
SumStatsSummary::accumulate(const vector<LocBin *> &loci)
{
    //
    // We are calculating the mean, variance, and standard deviation for several variables.
    //   We will calculate them partially, for each set of loci input to the program using
    //   the algorithm described in:
    //     B. P. Welford. (1962) Note on a Method for Calculating Corrected Sums of Squares and
    //     Products. Technometrics: 4(3), pp. 419-420.
    //
    CSLocus        *cloc;
    const LocSum   *s;
    const LocTally *t;

    for (uint i = 0; i < loci.size(); i++) {
        cloc = loci[i]->cloc;
        t    = loci[i]->s->meta_pop();

        size_t site_cnt = 0;

        for (uint pos = 0; pos < cloc->len; pos++) {
            //
            // Compile private alleles
            //
            if (t->nucs[pos].priv_allele >= 0)
                _private_cnt[t->nucs[pos].priv_allele]++;

            if (t->nucs[pos].allele_cnt == 2) {
                site_cnt++;

                for (uint pop = 0; pop < this->_pop_cnt; pop++) {

                    s = loci[i]->s->per_pop(pop);

                    if (s->nucs[pos].num_indv == 0) continue;

                    _n[pop]++;

                    if (s->nucs[pos].pi > 0) _var_sites[pop]++;

                    //
                    // Record if site deviates from HWE.
                    //
                    if (calc_hwp && s->nucs[pos].stat[2] < p_value_cutoff) _sig_hwe_dev[pop]++;

                    //
                    // Accumulate sums for each variable to calculate the means.
                    //
                    _num_indv_mean[pop] += s->nucs[pos].num_indv;
                    _p_mean[pop]        += s->nucs[pos].p;
                    _obs_het_mean[pop]  += s->nucs[pos].obs_het;
                    _obs_hom_mean[pop]  += s->nucs[pos].obs_hom;
                    _exp_het_mean[pop]  += s->nucs[pos].exp_het;
                    _exp_hom_mean[pop]  += s->nucs[pos].exp_hom;
                    _pi_mean[pop]       += s->nucs[pos].stat[0];
                    _fis_mean[pop]      += s->nucs[pos].stat[1] != -7.0 ? s->nucs[pos].stat[1] : 0.0;

                    _n_all[pop]++;
                    _num_indv_mean_all[pop] += s->nucs[pos].num_indv;
                    _p_mean_all[pop]        += s->nucs[pos].p;
                    _obs_het_mean_all[pop]  += s->nucs[pos].obs_het;
                    _obs_hom_mean_all[pop]  += s->nucs[pos].obs_hom;
                    _exp_het_mean_all[pop]  += s->nucs[pos].exp_het;
                    _exp_hom_mean_all[pop]  += s->nucs[pos].exp_hom;
                    _pi_mean_all[pop]       += s->nucs[pos].stat[0];
                    _fis_mean_all[pop]      += s->nucs[pos].stat[1] != -7.0 ? s->nucs[pos].stat[1] : 0.0;

                    //
                    // Accumulate a partial sum of squares to calculate the variance.
                    //
                    _num_indv_var[pop] += this->online_variance(s->nucs[pos].num_indv, _num_indv_acc_mean[pop], _n[pop]);
                    _p_var[pop]        += this->online_variance(s->nucs[pos].p,        _p_acc_mean[pop],        _n[pop]);
                    _obs_het_var[pop]  += this->online_variance(s->nucs[pos].obs_het,  _obs_het_acc_mean[pop],  _n[pop]);
                    _obs_hom_var[pop]  += this->online_variance(s->nucs[pos].obs_hom,  _obs_hom_acc_mean[pop],  _n[pop]);
                    _exp_het_var[pop]  += this->online_variance(s->nucs[pos].exp_het,  _exp_het_acc_mean[pop],  _n[pop]);
                    _exp_hom_var[pop]  += this->online_variance(s->nucs[pos].exp_hom,  _exp_hom_acc_mean[pop],  _n[pop]);
                    _pi_var[pop]       += this->online_variance(s->nucs[pos].stat[0],  _pi_acc_mean[pop],       _n[pop]);
                    _fis_var[pop]      += this->online_variance(s->nucs[pos].stat[1] != -7.0 ? s->nucs[pos].stat[1] : 0.0, _fis_acc_mean[pop], _n[pop]);

                    _num_indv_var_all[pop] += this->online_variance(s->nucs[pos].num_indv, _num_indv_acc_mean_all[pop], _n_all[pop]);
                    _p_var_all[pop]        += this->online_variance(s->nucs[pos].p,        _p_acc_mean_all[pop],        _n_all[pop]);
                    _obs_het_var_all[pop]  += this->online_variance(s->nucs[pos].obs_het,  _obs_het_acc_mean_all[pop],  _n_all[pop]);
                    _obs_hom_var_all[pop]  += this->online_variance(s->nucs[pos].obs_hom,  _obs_hom_acc_mean_all[pop],  _n_all[pop]);
                    _exp_het_var_all[pop]  += this->online_variance(s->nucs[pos].exp_het,  _exp_het_acc_mean_all[pop],  _n_all[pop]);
                    _exp_hom_var_all[pop]  += this->online_variance(s->nucs[pos].exp_hom,  _exp_hom_acc_mean_all[pop],  _n_all[pop]);
                    _pi_var_all[pop]       += this->online_variance(s->nucs[pos].stat[0],  _pi_acc_mean_all[pop],       _n_all[pop]);
                    _fis_var_all[pop]      += this->online_variance(s->nucs[pos].stat[1] != -7.0 ? s->nucs[pos].stat[1] : 0.0, _fis_acc_mean_all[pop], _n_all[pop]);
                }

            } else if (t->nucs[pos].allele_cnt == 1) {
                site_cnt++;

                for (uint pop = 0; pop < this->_pop_cnt; pop++) {
                    s = loci[i]->s->per_pop(pop);

                    if (s->nucs[pos].num_indv == 0) continue;

                    _n_all[pop]++;
                    _num_indv_mean_all[pop] += s->nucs[pos].num_indv;
                    _p_mean_all[pop]        += s->nucs[pos].p;
                    _obs_het_mean_all[pop]  += s->nucs[pos].obs_het;
                    _obs_hom_mean_all[pop]  += s->nucs[pos].obs_hom;
                    _exp_het_mean_all[pop]  += s->nucs[pos].exp_het;
                    _exp_hom_mean_all[pop]  += s->nucs[pos].exp_hom;
                    _pi_mean_all[pop]       += s->nucs[pos].stat[0];
                    _fis_mean_all[pop]      += s->nucs[pos].stat[1] != -7.0 ? s->nucs[pos].stat[1] : 0.0;

                    _num_indv_var_all[pop] += this->online_variance(s->nucs[pos].num_indv, _num_indv_acc_mean_all[pop], _n_all[pop]);
                    _p_var_all[pop]        += this->online_variance(s->nucs[pos].p,        _p_acc_mean_all[pop],        _n_all[pop]);
                    _obs_het_var_all[pop]  += this->online_variance(s->nucs[pos].obs_het,  _obs_het_acc_mean_all[pop],  _n_all[pop]);
                    _obs_hom_var_all[pop]  += this->online_variance(s->nucs[pos].obs_hom,  _obs_hom_acc_mean_all[pop],  _n_all[pop]);
                    _exp_het_var_all[pop]  += this->online_variance(s->nucs[pos].exp_het,  _exp_het_acc_mean_all[pop],  _n_all[pop]);
                    _exp_hom_var_all[pop]  += this->online_variance(s->nucs[pos].exp_hom,  _exp_hom_acc_mean_all[pop],  _n_all[pop]);
                    _pi_var_all[pop]       += this->online_variance(s->nucs[pos].stat[0],  _pi_acc_mean_all[pop],       _n_all[pop]);
                    _fis_var_all[pop]      += this->online_variance(s->nucs[pos].stat[1] != -7.0 ? s->nucs[pos].stat[1] : 0.0, _fis_acc_mean_all[pop], _n_all[pop]);
                }
            }
        }

        //
        // Accumulate the locus length and overlap.
        //
        _locus_n++;
        _locus_gt_sites_mean += site_cnt;
        _locus_gt_sites_var  += this->online_variance(site_cnt, _locus_gt_sites_acc_mean, _locus_n);

        if (cloc->pe_ctg) {
            _locus_pe_ctg_n++;
            _locus_len_mean_all += cloc->len - 10;
            _locus_len_var_all  += this->online_variance(cloc->len - 10, _locus_len_acc_mean_all, _locus_pe_ctg_n);

            if (cloc->overlap > 0) {
                _locus_overlap_n++;
                _locus_len_mean += cloc->len;
                _locus_len_var  += this->online_variance(cloc->len, _locus_len_acc_mean, _locus_overlap_n);
                _overlap_mean   += cloc->overlap;
                _overlap_var    += this->online_variance(cloc->overlap, _overlap_acc_mean, _locus_overlap_n);
            }
        }
    }

    return 0;
}

inline double
SumStatsSummary::online_variance(double x, double &acc_mean, double n)
{
    double delta1, delta2;

    delta1    = x - acc_mean;
    acc_mean += delta1 / n;
    delta2    = x - acc_mean;
    return delta1 * delta2;
}

int
SumStatsSummary::final_calculation()
{
    //
    // Finish the mean calculation.
    //
    for (uint j = 0; j < this->_pop_cnt; j++) {
        _num_indv_mean[j] = _num_indv_mean[j] / _n[j];
        _p_mean[j]        = _p_mean[j]        / _n[j];
        _obs_het_mean[j]  = _obs_het_mean[j]  / _n[j];
        _obs_hom_mean[j]  = _obs_hom_mean[j]  / _n[j];
        _exp_het_mean[j]  = _exp_het_mean[j]  / _n[j];
        _exp_hom_mean[j]  = _exp_hom_mean[j]  / _n[j];
        _pi_mean[j]       = _pi_mean[j]       / _n[j];
        _fis_mean[j]      = _fis_mean[j]      / _n[j];

        _num_indv_mean_all[j] = _num_indv_mean_all[j] / _n_all[j];
        _p_mean_all[j]        = _p_mean_all[j]        / _n_all[j];
        _obs_het_mean_all[j]  = _obs_het_mean_all[j]  / _n_all[j];
        _obs_hom_mean_all[j]  = _obs_hom_mean_all[j]  / _n_all[j];
        _exp_het_mean_all[j]  = _exp_het_mean_all[j]  / _n_all[j];
        _exp_hom_mean_all[j]  = _exp_hom_mean_all[j]  / _n_all[j];
        _pi_mean_all[j]       = _pi_mean_all[j]       / _n_all[j];
        _fis_mean_all[j]      = _fis_mean_all[j]      / _n_all[j];
    }

    _locus_len_mean      = _locus_len_mean      / _locus_overlap_n;
    _locus_len_mean_all  = _locus_len_mean_all  / _locus_pe_ctg_n;
    _overlap_mean        = _overlap_mean        / _locus_overlap_n;
    _locus_gt_sites_mean = _locus_gt_sites_mean / _locus_n;

    //
    // Finish the online variance calculation.
    //
    for (uint j = 0; j < this->_pop_cnt; j++) {
        _num_indv_var[j] = _num_indv_var[j] / (_n[j] - 1);
        _p_var[j]        = _p_var[j]        / (_n[j] - 1);
        _obs_het_var[j]  = _obs_het_var[j]  / (_n[j] - 1);
        _obs_hom_var[j]  = _obs_hom_var[j]  / (_n[j] - 1);
        _exp_het_var[j]  = _exp_het_var[j]  / (_n[j] - 1);
        _exp_hom_var[j]  = _exp_hom_var[j]  / (_n[j] - 1);
        _pi_var[j]       = _pi_var[j]       / (_n[j] - 1);
        _fis_var[j]      = _fis_var[j]      / (_n[j] - 1);

        _num_indv_var_all[j] = _num_indv_var_all[j] / (_n_all[j] - 1);
        _p_var_all[j]        = _p_var_all[j]        / (_n_all[j] - 1);
        _obs_het_var_all[j]  = _obs_het_var_all[j]  / (_n_all[j] - 1);
        _obs_hom_var_all[j]  = _obs_hom_var_all[j]  / (_n_all[j] - 1);
        _exp_het_var_all[j]  = _exp_het_var_all[j]  / (_n_all[j] - 1);
        _exp_hom_var_all[j]  = _exp_hom_var_all[j]  / (_n_all[j] - 1);
        _pi_var_all[j]       = _pi_var_all[j]       / (_n_all[j] - 1);
        _fis_var_all[j]      = _fis_var_all[j]      / (_n_all[j] - 1);
    }

    _locus_len_var     = _locus_len_var     / (_locus_overlap_n - 1);
    _locus_len_var_all = _locus_len_var_all / (_locus_pe_ctg_n  - 1);
    _overlap_var       = _overlap_var       / (_locus_overlap_n - 1);

    _locus_gt_sites_var = _locus_gt_sites_var / (_locus_n - 1);

    //
    // Calculate the first half of the standard deviation.
    //
    for (uint j = 0; j < this->_pop_cnt; j++) {
        _sq_n[j]     = sqrt(_n[j]);
        _sq_n_all[j] = sqrt(_n_all[j]);
    }

    return 0;
}

int
SumStatsSummary::write_results()
{
    string   path = out_path + out_prefix + ".sumstats_summary.tsv";
    ofstream fh(path.c_str(), ofstream::out);
    if (fh.fail()) {
        cerr << "Error opening sumstats summary file '" << path << "'\n";
        exit(1);
    }
    fh.precision(fieldw);
    fh.setf(std::ios::fixed);

    //
    // Write out locus length and overlap statistics for de novo data.
    //
    ostream os (cout.rdbuf());
    os << std::fixed << std::setprecision(2);
    if (!loci_ordered && this->_locus_pe_ctg_n > 0) {
        os   << "Number of loci with PE contig: " << this->_locus_pe_ctg_n << " ("
             << as_percentage(this->_locus_pe_ctg_n, this->_locus_n) << ");\n"
             << "  Mean length of loci: " << this->_locus_len_mean_all << "bp "
             << "(stderr " << sqrt(this->_locus_len_var_all) / sqrt(this->_locus_n) << ");\n"
             << "Number of loci with SE/PE overlap: " << this->_locus_overlap_n << " ("
             << as_percentage(this->_locus_overlap_n, this->_locus_n) << ");\n"
             << "  Mean length of overlapping loci: " << this->_locus_len_mean << "bp "
             << "(stderr " << sqrt(this->_locus_len_var) / sqrt(this->_locus_n) << "); "
             << "mean overlap: " << this->_overlap_mean
             << "bp (stderr " << sqrt(this->_overlap_var) / sqrt(this->_locus_n) << ");\n";
    }
    os   << "Mean genotyped sites per locus: " << this->_locus_gt_sites_mean << "bp "
         << "(stderr " << sqrt(this->_locus_gt_sites_var) / sqrt(this->_locus_n) << ").\n";

    //
    // Write out summary statistics of the summary statistics.
    //
    fh << "# Variant positions\n"
       << "# Pop ID\t"
       << "Private\t"
       << "Num_Indv\t"
       << "Var\t"
       << "StdErr\t"
       << "P\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs_Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs_Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp_Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp_Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Pi\t"
       << "Var\t"
       << "StdErr\t"
       << "Fis\t"
       << "Var\t"
       << "StdErr\n";

    for (uint j = 0; j < this->_pop_cnt; j++)
        fh << mpopi.pops()[j].name << "\t"
           << _private_cnt[j]         << "\t"
           << _num_indv_mean[j]       << "\t"
           << _num_indv_var[j]        << "\t"
           << sqrt(_num_indv_var[j]) / _sq_n[j] << "\t"
           << _p_mean[j]              << "\t"
           << _p_var[j]               << "\t"
           << sqrt(_p_var[j])         / _sq_n[j] << "\t"
           << _obs_het_mean[j]        << "\t"
           << _obs_het_var[j]         << "\t"
           << sqrt(_obs_het_var[j])   / _sq_n[j] << "\t"
           << _obs_hom_mean[j]        << "\t"
           << _obs_hom_var[j]         << "\t"
           << sqrt(_obs_hom_var[j])   / _sq_n[j] << "\t"
           << _exp_het_mean[j]        << "\t"
           << _exp_het_var[j]         << "\t"
           << sqrt(_exp_het_var[j])   / _sq_n[j] << "\t"
           << _exp_hom_mean[j]        << "\t"
           << _exp_hom_var[j]         << "\t"
           << sqrt(_exp_hom_var[j])   / _sq_n[j] << "\t"
           << _pi_mean[j]             << "\t"
           << _pi_var[j]              << "\t"
           << sqrt(_pi_var[j])        / _sq_n[j] << "\t"
           << _fis_mean[j]            << "\t"
           << _fis_var[j]             << "\t"
           << sqrt(_num_indv_var[j])  / _sq_n[j] << "\n";

    cout << "\nPopulation summary statistics (more detail in populations.sumstats_summary.tsv):\n";

    for (uint j = 0; j < this->_pop_cnt; j++)
        cout << "  " << mpopi.pops()[j].name << ": "
             << setprecision(fieldw) << _num_indv_mean[j] << " samples per locus; "
             << "pi: " << _pi_mean[j] << "; "
             << setprecision(10) << "all/variant/polymorphic sites: " << _n_all[j] << "/" << _n[j] << "/" << _var_sites[j] << "; "
             << setprecision(fieldw) << "private alleles: " << _private_cnt[j] << "\n";

    fh << "# All positions (variant and fixed)\n"
       << "# Pop ID\t"
       << "Private\t"
       << "Sites\t"
       << "Variant_Sites\t"
       << "Polymorphic_Sites\t"
       << "%Polymorphic_Loci\t"
       << "Num_Indv\t"
       << "Var\t"
       << "StdErr\t"
       << "P\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs_Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs_Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp_Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp_Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Pi\t"
       << "Var\t"
       << "StdErr\t"
       << "Fis\t"
       << "Var\t"
       << "StdErr\n";

    for (uint j = 0; j < this->_pop_cnt; j++) {
        fh << mpopi.pops()[j].name << "\t"
           << _private_cnt[j]             << "\t" << setprecision(0)
           << _n_all[j]                   << "\t"
           << _n[j]                       << "\t"
           << _var_sites[j]               << "\t" << setprecision(fieldw)
           << _var_sites[j]             / _n_all[j] * 100 << "\t"
           << _num_indv_mean_all[j]       << "\t"
           << _num_indv_var_all[j]        << "\t"
           << sqrt(_num_indv_var_all[j]) / _sq_n_all[j] << "\t"
           << _p_mean_all[j]              << "\t"
           << _p_var_all[j]               << "\t"
           << sqrt(_p_var_all[j])        / _sq_n_all[j] << "\t"
           << _obs_het_mean_all[j]        << "\t"
           << _obs_het_var_all[j]         << "\t"
           << sqrt(_obs_het_var_all[j])  / _sq_n_all[j] << "\t"
           << _obs_hom_mean_all[j]        << "\t"
           << _obs_hom_var_all[j]         << "\t"
           << sqrt(_obs_hom_var_all[j])  / _sq_n_all[j] << "\t"
           << _exp_het_mean_all[j]        << "\t"
           << _exp_het_var_all[j]         << "\t"
           << sqrt(_exp_het_var_all[j])  / _sq_n_all[j] << "\t"
           << _exp_hom_mean_all[j]        << "\t"
           << _exp_hom_var_all[j]         << "\t"
           << sqrt(_exp_hom_var_all[j])  / _sq_n_all[j] << "\t"
           << _pi_mean_all[j]             << "\t"
           << _pi_var_all[j]              << "\t"
           << sqrt(_pi_var_all[j])       / _sq_n_all[j] << "\t"
           << _fis_mean_all[j]            << "\t"
           << _fis_var_all[j]             << "\t"
           << sqrt(_num_indv_var_all[j]) / _sq_n_all[j] << "\n";
    }

    fh.close();

    if (calc_hwp) {
        os << "\n"
           << "Number of variable sites found to be significantly out of Hardy-Weinberg equilibrium (<" << p_value_cutoff << "):\n";
        for (uint j = 0; j < this->_pop_cnt; j++)
            os << "  "
               << mpopi.pops()[j].name << ": "
               << _sig_hwe_dev[j]      << "\n";
    }

    return 0;
}

SumStatsSummary::~SumStatsSummary()
{
    delete [] this->_private_cnt;
    delete [] this->_n;
    delete [] this->_var_sites;
    delete [] this->_sq_n;
    delete [] this->_num_indv_mean;
    delete [] this->_num_indv_acc_mean;
    delete [] this->_num_indv_var;
    delete [] this->_p_mean;
    delete [] this->_p_acc_mean;
    delete [] this->_p_var;
    delete [] this->_obs_het_mean;
    delete [] this->_obs_het_acc_mean;
    delete [] this->_obs_het_var;
    delete [] this->_obs_hom_mean;
    delete [] this->_obs_hom_acc_mean;
    delete [] this->_obs_hom_var;
    delete [] this->_exp_het_mean;
    delete [] this->_exp_het_acc_mean;
    delete [] this->_exp_het_var;
    delete [] this->_exp_hom_mean;
    delete [] this->_exp_hom_acc_mean;
    delete [] this->_exp_hom_var;
    delete [] this->_pi_mean;
    delete [] this->_pi_acc_mean;
    delete [] this->_pi_var;
    delete [] this->_fis_mean;
    delete [] this->_fis_acc_mean;
    delete [] this->_fis_var;

    delete [] this->_n_all;
    delete [] this->_sq_n_all;
    delete [] this->_num_indv_mean_all;
    delete [] this->_num_indv_acc_mean_all;
    delete [] this->_num_indv_var_all;
    delete [] this->_p_mean_all;
    delete [] this->_p_acc_mean_all;
    delete [] this->_p_var_all;
    delete [] this->_obs_het_mean_all;
    delete [] this->_obs_het_acc_mean_all;
    delete [] this->_obs_het_var_all;
    delete [] this->_obs_hom_mean_all;
    delete [] this->_obs_hom_acc_mean_all;
    delete [] this->_obs_hom_var_all;
    delete [] this->_exp_het_mean_all;
    delete [] this->_exp_het_acc_mean_all;
    delete [] this->_exp_het_var_all;
    delete [] this->_exp_hom_mean_all;
    delete [] this->_exp_hom_acc_mean_all;
    delete [] this->_exp_hom_var_all;
    delete [] this->_pi_mean_all;
    delete [] this->_pi_acc_mean_all;
    delete [] this->_pi_var_all;
    delete [] this->_fis_mean_all;
    delete [] this->_fis_acc_mean_all;
    delete [] this->_fis_var_all;
}

int
correct_fst_bonferroni_win(vector<PopPair *> &pairs)
{
    int      limit = 3 * sigma;
    int      limit_l, limit_u;
    double   correction;
    uint     cnt, pos_l, pos_u;

    pos_l = 0;
    pos_u = 0;

    for (uint pos_c = 0; pos_c < pairs.size(); pos_c++) {
        if (pairs[pos_c] == NULL) continue;

        limit_l = pairs[pos_c]->bp - limit > 0 ? pairs[pos_c]->bp - limit : 0;
        limit_u = pairs[pos_c]->bp + limit;

        while (pos_l <  pairs.size()) {
            if (pairs[pos_l] == NULL) {
                pos_l++;
            } else {
                if (pairs[pos_l]->bp < limit_l)
                    pos_l++;
                else
                    break;
            }
        }
        while (pos_u < pairs.size()) {
            if (pairs[pos_u] == NULL) {
                pos_u++;
            } else {
                if (pairs[pos_u]->bp < limit_u)
                    pos_u++;
                else
                    break;
            }
        }

        cnt = 0;
        for (uint i = pos_l; i < pos_u; i++) {
            if (pairs[i] == NULL) continue;
            cnt++;
        }

        correction = p_value_cutoff / cnt;
        pairs[pos_c]->stat[0] = pairs[pos_c]->fet_p < correction ? pairs[pos_c]->fst : 0;
    }

    return 0;
}

int
LocusSmoothing::snpstats(const vector<LocBin *> &loci, ofstream &log_fh)
{
    for (uint i = 0; i < this->_mpopi->pops().size(); i++) {

        vector<const SumStat *> sites;

        this->_ord_ss->order(sites, loci, i);
        this->_ks_ss->smooth(sites);
    }

    return 0;
}

int
LocusSmoothing::hapstats(const vector<LocBin *> &loci, ofstream &log_fh)
{
    vector<const LocStat *> sites;

    for (uint i = 0; i < this->_mpopi->pops().size(); i++) {
        sites.clear();

        this->_ord_ls->order(sites, loci);
        this->_ks_ls->smooth(sites);
    }

    return 0;
}

int
LocusSmoothing::snp_divergence(const vector<LocBin *> &loci, const vector<vector<PopPair **>> &div, ofstream &log_fh)
{
    vector<const PopPair *> sites;

    for (uint i = 0; i < div.size(); i++) {
        assert(div[i].size() == loci.size());

        sites.clear();

        if (this->_ord_pp->order(sites, loci, div[i]))
            this->_ks_pp->smooth(sites);
    }

    return 0;
}

int
LocusSmoothing::hap_divergence(const vector<LocBin *> &loci,
                               const vector<vector<HapStat *>> &div,
                               const vector<HapStat *> &metadiv,
                               ofstream &log_fh)
{
    vector<const HapStat *> sites;

    for (uint i = 0; i < div.size(); i++) {
        assert(div[i].size() == loci.size());

        sites.clear();

        this->_ord_hs->order(sites, loci, div[i]);
        this->_ks_hs->smooth(sites);
    }

    //
    // Kernel-smooth the haplotype divergence statistics for the metapopulation.
    //
    sites.clear();

    this->_ord_hs->order(sites, loci, metadiv);
    this->_ks_hs->smooth(sites);

    return 0;
}

int
bootstrap_popstats_approximate_dist(vector<double> &fis_samples,
                                    vector<double> &pi_samples,
                                    vector<int>  &allele_samples,
                                    double *weights, int *snp_dist, int sites_per_snp,
                                    map<int, vector<double> > &approx_fis_dist,
                                    map<int, vector<double> > &approx_pi_dist)
{
    //
    // Allocate an array of bootstrap resampling objects.
    //
    int win_size = 6 * sigma + 1;
    int win_cntr = win_size / 2;

    //
    // Initialize the Fst distribution map.
    //
    for (int i = 0; i < max_snp_dist; i++) {
        if (snp_dist[i] == 0.0) continue;

        // cout << "SNP Dist: " << i << " snps occurred " << snp_dist[i] << "\n";
        approx_fis_dist[i] = vector<double> ();
        approx_fis_dist[i].reserve(bootstrap_reps);

        approx_pi_dist[i] = vector<double> ();
        approx_pi_dist[i].reserve(bootstrap_reps);
    }

    vector<int> poss;
    poss.reserve(max_snp_dist);
    double weighted_fis, weighted_pi, sum_fis, sum_pi, final_weight_fis, final_weight_pi;
    // int    index_1, index_2;
    int    pos, index_3, dist, start, end;
    int    half = sites_per_snp / 2;

    for (int i = 0; i < max_snp_dist; i++) {
        if (snp_dist[i] == 0.0) continue;

        cout << "  Generating NULL distribution for " << i << " SNPs...\n";

        // #pragma omp parallel private(poss, pos, index_1, index_2, index_3, dist, sum_fis, sum_pi, weighted_fis, weighted_pi, final_weight_fis, final_weight_pi)
        #pragma omp parallel private(poss, pos, index_3, dist, sum_fis, sum_pi, weighted_fis, weighted_pi, final_weight_fis, final_weight_pi)
        {
            BSample *bs  = new BSample[win_size];

            //
            // Populate the BSample objects.
            //
            for (int n = 0; n < win_size;  n++)
                bs[n].bp = n + 1;

            vector<double> fiss, pis;

            //
            // Bootstrap this bitch.
            //
            #pragma omp for schedule(dynamic, 1)
            for (int j = 0; j < bootstrap_reps; j++) {
                // cout << "    Bootsrap rep " << j << "\n";

                //
                // First SNP is always placed at the center of the window.
                //
                pos     = win_cntr;
                // index_1 = (int) (fis_samples.size()    * (random() / (RAND_MAX + 1.0)));
                // index_2 = (int) (pi_samples.size()     * (random() / (RAND_MAX + 1.0)));
                index_3 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));
                //
                // Fill in the area around the SNP with fixed sites.
                //
                start = pos - half > 0 ? pos - half : 0;
                end   = pos + half < win_size ? pos + half : win_size;
                for (int n = start; n < end; n++) {
                    // bs[n].f       = 0;
                    // bs[n].pi      = 0;
                    bs[n].alleles = bs[pos].alleles;
                    poss.push_back(n);
                }
                // bs[pos].f       = fis_samples[index_1];
                // bs[pos].pi      = pi_samples[index_2];
                bs[pos].alleles = allele_samples[index_3];
                // cout << "      Placing SNP at position: " << pos << "; with data from " << index_1 << " filling area from " << start << " to " << end << "\n";

                //
                // Randomly select the positions and values for each SNP to populate the window
                //
                for (int k = 0; k < i - 1; k++) {
                    pos     = (int) (win_size * (random() / (RAND_MAX + 1.0)));
                    // index_1 = (int) (fis_samples.size()    * (random() / (RAND_MAX + 1.0)));
                    // index_2 = (int) (pi_samples.size()     * (random() / (RAND_MAX + 1.0)));
                    index_3 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));

                    poss.push_back(pos);
                    //
                    // Fill in the area around the SNP with fixed sites.
                    //
                    start = pos - half > 0 ? pos - half : 0;
                    end   = pos + half < win_size ? pos + half : win_size;
                    for (int n = start; n < end; n++) {
                        // bs[n].f       = 0;
                        // bs[n].pi      = 0;
                        bs[n].alleles = bs[pos].alleles;
                        poss.push_back(n);
                    }
                    // bs[pos].f       = fis_samples[index_1];
                    // bs[pos].pi      = pi_samples[index_2];
                    bs[pos].alleles = allele_samples[index_3];
                    // cout << "      Placing SNP at position: " << pos << "; with data from " << index_1 << " filling area from " << start << " to " << end << "\n";
                }

                weighted_fis = 0.0;
                sum_fis      = 0.0;
                weighted_pi  = 0.0;
                sum_pi       = 0.0;

                for (int n = 0; n < win_size; n++) {
                    // if (bs[n].pi < 0.0)
                    // continue;
                    //
                    // Calculate weighted Fst at this position.
                    //
                    dist = bs[n].bp > bs[win_cntr].bp ? bs[n].bp - bs[win_cntr].bp : bs[win_cntr].bp - bs[n].bp;

                    final_weight_fis = (bs[n].alleles - 1) * weights[dist];
                    // weighted_fis    += bs[n].f * final_weight_fis;
                    sum_fis         += final_weight_fis;

                    final_weight_pi  = (bs[n].alleles - 1) * weights[dist];
                    // weighted_pi     += bs[n].pi * final_weight_pi;
                    sum_pi          += final_weight_pi;
                }

                fiss.push_back(weighted_fis / sum_fis);
                pis.push_back(weighted_pi  / sum_pi);
                // cout << "      New weighted fis value: " << weighted_fis / sum_fis << "; size: " << fiss.size() << "\n";

                for (uint n = 0; n < poss.size(); n++) {
                    // bs[poss[n]].f  = 0.0;
                    // bs[poss[n]].pi = -1.0;
                }
                poss.clear();
            }

//          #pragma omp critical
//          {
//              vector<double> &f = approx_fis_dist[i];
//              for (uint n = 0; n < fiss.size(); n++)
//                  f.push_back(fiss[n]);
//              vector<double> &p = approx_pi_dist[i];
//              for (uint n = 0; n < pis.size(); n++)
//                  p.push_back(pis[n]);
//          }

            delete [] bs;
        }

        sort(approx_fis_dist[i].begin(), approx_fis_dist[i].end());
        sort(approx_pi_dist[i].begin(),  approx_pi_dist[i].end());
    }

    return 0;
}

int
bootstrap_fst_approximate_dist(vector<double> &fst_samples,
                               vector<int>  &allele_samples,
                               double *weights, int *snp_dist,
                               map<int, vector<double> > &approx_fst_dist)
{
    //
    // Allocate an array of bootstrap resampling objects.
    //
    int win_size = 6 * sigma + 1;
    int win_cntr = win_size / 2;

    //
    // Initialize the Fst distribution map.
    //
    for (int i = 0; i < max_snp_dist; i++) {
        if (snp_dist[i] == 0.0) continue;

        // cout << "SNP Dist: " << i << " snps occurred " << snp_dist[i] << "\n";
        approx_fst_dist[i] = vector<double> ();
        approx_fst_dist[i].reserve(bootstrap_reps);
    }

    vector<int> poss;
    poss.reserve(max_snp_dist);
    double weighted_fst, sum, final_weight;
    //int    index_1;
    int    pos, index_2, dist;

    for (int i = 0; i < max_snp_dist; i++) {
        if (snp_dist[i] == 0.0) continue;

        cout << "  Generating NULL distribution for " << i << " SNPs...\n";

        // #pragma omp parallel private(poss, pos, index_1, index_2, dist, sum, weighted_fst, final_weight)
        #pragma omp parallel private(poss, pos, index_2, dist, sum, weighted_fst, final_weight)
        {
            BSample *bs  = new BSample[win_size];

            //
            // Populate the BSample objects.
            //
            for (int n = 0; n < win_size;  n++)
                bs[n].bp = n + 1;

            vector<double> fsts;

            //
            // Bootstrap this bitch.
            //
            #pragma omp for schedule(dynamic, 1)
            for (int j = 0; j < bootstrap_reps; j++) {
                // cout << "Bootsrap rep " << j << "\n";

                //
                // First SNP is always placed at the center of the window.
                //
                pos     = win_cntr;
                // index_1 = (int) (fst_samples.size() * (random() / (RAND_MAX + 1.0)));
                index_2 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));
                // bs[pos].f       = fst_samples[index_1];
                bs[pos].alleles = allele_samples[index_2];

                //
                // Randomly select the positions and values for each SNP to populate the window
                //
                for (int k = 0; k < i - 1; k++) {
                    pos     = (int) (win_size * (random() / (RAND_MAX + 1.0)));
                    // index_1 = (int) (fst_samples.size() * (random() / (RAND_MAX + 1.0)));
                    index_2 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));
                    // bs[pos].f       = fst_samples[index_1];
                    // bs[pos].alleles = allele_samples[index_2];
                    // cout << "  " << j << ": Placing SNP at position: " << pos << " with data from index " << index_1 << "\n";

                    poss.push_back(pos);
                }

                weighted_fst = 0.0;
                sum          = 0.0;

                for (int n = 0; n < win_size; n++) {
                    // if (bs[n].f == 0.0)
                    // continue;
                    //
                    // Calculate weighted Fst at this position.
                    //
                    dist = bs[n].bp > bs[win_cntr].bp ? bs[n].bp - bs[win_cntr].bp : bs[win_cntr].bp - bs[n].bp;

                    final_weight  = (bs[n].alleles - 1) * weights[dist];
                    // weighted_fst += bs[n].f * final_weight;
                    sum          += final_weight;
                }

                fsts.push_back(weighted_fst / sum);
                // cout << "    New weighted Fst value: " << weighted_fst / sum << "; size: " << fsts.size() << "\n";

                // for (uint n = 0; n < poss.size(); n++)
                // bs[poss[n]].f = 0.0;
                poss.clear();
            }

//          #pragma omp critical
//          {
//              vector<double> &f = approx_fst_dist[i];
//              for (uint n = 0; n < fsts.size(); n++)
//                  f.push_back(fsts[n]);
//          }

            delete [] bs;
        }

        sort(approx_fst_dist[i].begin(), approx_fst_dist[i].end());
    }

    return 0;
}

double
bootstrap_approximate_pval(int snp_cnt, double stat, map<int, vector<double> > &approx_dist)
{
    if (approx_dist.count(snp_cnt) == 0)
        return 1.0;

    vector<double>::iterator up;
    vector<double> &dist = approx_dist[snp_cnt];
    double pos;

    up  = upper_bound(dist.begin(), dist.end(), stat);

    if (up == dist.begin())
        pos = 1;
    else if (up == dist.end())
        pos = dist.size();
    else
        pos = up - dist.begin() + 1;

    double res = 1.0 - (pos / (double) dist.size());

    // cout << "Generated Approx Smoothed Fst Distribution:\n";
    // for (uint n = 0; n < dist.size(); n++)
    //  cout << "  n: " << n << "; Fst: " << dist[n] << "\n";

    // cout << "Comparing Fst value: " << stat
    //   << " at position " << (up - dist.begin()) << " out of "
    //   << dist.size() << " positions (converted position: " << pos << "); pvalue: " << res << ".\n";

    return res;
}

int
LocusFilter::load_blacklist(string path)
{
    char     line[id_len];
    ifstream fh(path.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening white/black list file '" << path << "'\n";
        exit(1);
    }

    size_t line_num = 0;
    while (fh.getline(line, id_len)) {
        ++line_num;

        //
        // Skip blank & commented lines ; correct windows-style line ends.
        //
        size_t len = strlen(line);
        if (len == 0) {
            continue;
        } else if (line[len-1] == '\r') {
            line[len-1] = '\0';
            --len;
            if (len == 0)
                continue;
        }
        char* p = line;
        while (isspace(*p) && *p != '\0')
            ++p;
        if (*p == '#')
            continue;

        //
        // Parse the blacklist
        //
        char* e;
        int marker = (int) strtol(line, &e, 10);
        if (*e == '\0') {
            this->_blacklist.insert(marker);
        } else {
            cerr << "Error: Unable to parse blacklist '" << path << "' at line " << line_num << ".\n";
            throw exception();
        }
    }

    fh.close();

    if (this->_blacklist.size() == 0) {
        cerr << "Error: Unable to load any markers from '" << path << "'\n";
        exit(1);
    }

    return (int) this->_blacklist.size();
}

int load_marker_list(string path, set<int> &list) {
    char     line[id_len];
    ifstream fh(path.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening white/black list file '" << path << "'\n";
        exit(1);
    }

    size_t line_num = 0;
    while (fh.getline(line, id_len)) {
        ++line_num;

        //
        // Skip blank & commented lines ; correct windows-style line ends.
        //
        size_t len = strlen(line);
        if (len == 0) {
            continue;
        } else if (line[len-1] == '\r') {
            line[len-1] = '\0';
            --len;
            if (len == 0)
                continue;
        }
        char* p = line;
        while (isspace(*p) && *p != '\0')
            ++p;
        if (*p == '#')
            continue;

        //
        // Parse the blacklist
        //
        char* e;
        int marker = (int) strtol(line, &e, 10);
        if (*e == '\0') {
            list.insert(marker);
        } else {
            cerr << "Error: Unable to parse blacklist '" << path << "' at line " << line_num << ".\n";
            throw exception();
        }
    }

    fh.close();

    if (list.size() == 0) {
        cerr << "Error: Unable to load any markers from '" << path << "'\n";
        exit(1);
    }

    return 0;
}

int
LocusFilter::load_whitelist(string path)
{
    char     line[id_len];
    ifstream fh(path.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening white/black list file '" << path << "'\n";
        exit(1);
    }

    vector<string> parts;
    uint col;
    char *e;

    uint line_num = 0;
    while (fh.getline(line, id_len)) {
        ++line_num;

        //
        // Skip blank & commented lines ; correct windows-style line ends.
        //
        size_t len = strlen(line);
        if (len == 0) {
            continue;
        } else if (line[len-1] == '\r') {
            line[len-1] = '\0';
            --len;
            if (len == 0)
                continue;
        }
        char* p = line;
        while (isspace(*p) && *p != '\0')
            ++p;
        if (*p == '#')
            continue;

        //
        // Parse the whitelist, we expect:
        // <marker>[<tab><snp column>]
        //
        parse_tsv(line, parts);

        if (parts.size() > 2) {
            cerr << "Error: Too many columns in whitelist " << path << "' at line " << line_num << "\n";
            exit(1);
        }
        int marker = (int) strtol(parts[0].c_str(), &e, 10);
        if (*e != '\0') {
            cerr << "Error: Unable to parse whitelist, '" << path << "' at line " << line_num << "\n";
            exit(1);
        }
        this->_whitelist[marker];
        if (parts.size() == 2) {
            col = (int) strtol(parts[1].c_str(), &e, 10);
            if (*e != '\0') {
                cerr << "Error: Unable to parse whitelist, '" << path << "' at line " << line_num << "\n";
                exit(1);
            }
            this->_whitelist[marker].insert(col);
        }
    }
    if (this->_whitelist.size() == 0) {
        cerr << "Error: Unable to load any markers from '" << path << "'\n";
        help();
    }
    return (int) this->_whitelist.size();
}

int load_marker_column_list(string path, map<int, set<int> > &list) {
    char     line[id_len];
    ifstream fh(path.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening white/black list file '" << path << "'\n";
        exit(1);
    }

    vector<string> parts;
    uint col;
    char *e;

    uint line_num = 1;
    while (fh.getline(line, id_len)) {

        //
        // Skip blank & commented lines ; correct windows-style line ends.
        //
        size_t len = strlen(line);
        if (len == 0) {
            continue;
        } else if (line[len-1] == '\r') {
            line[len-1] = '\0';
            --len;
            if (len == 0)
                continue;
        }
        char* p = line;
        while (isspace(*p) && *p != '\0')
            ++p;
        if (*p == '#')
            continue;

        //
        // Parse the whitelist, we expect:
        // <marker>[<tab><snp column>]
        //
        parse_tsv(line, parts);

        if (parts.size() > 2) {
            cerr << "Error: Too many columns in whitelist " << path << "' at line " << line_num << "\n";
            exit(1);

        } else if (parts.size() == 2) {
            int marker = (int) strtol(parts[0].c_str(), &e, 10);
            if (*e != '\0') {
                cerr << "Error: Unable to parse whitelist, '" << path << "' at line " << line_num << "\n";
                exit(1);
            }
            col = (int) strtol(parts[1].c_str(), &e, 10);
            if (*e != '\0') {
                cerr << "Error: Unable to parse whitelist, '" << path << "' at line " << line_num << "\n";
                exit(1);
            }
            list[marker].insert(col);

        } else {
            int marker = (int) strtol(parts[0].c_str(), &e, 10);
            if (*e != '\0') {
                cerr << "Error: Unable to parse whitelist, '" << path << "' at line " << line_num << "\n";
                exit(1);
            }
            list.insert(make_pair(marker, set<int>()));
        }

        line_num++;
    }

    fh.close();

    if (list.size() == 0) {
        cerr << "Error: Unable to load any markers from '" << path << "'\n";
        help();
    }

    return 0;
}

bool
hap_compare(const pair<string, int>& a, const pair<string, int>& b)
{
    return (a.second > b.second);
}

void
output_parameters(ostream &fh)
{
    fh << "populations parameters selected:\n";
    if (input_mode == InputMode::vcf)
        fh << "  Input mode: VCF\n";
    fh
        << "  Percent samples limit per population: " << min_samples_per_pop << "\n"
        << "  Locus Population limit: " << min_populations << "\n"
        << "  Percent samples overall: " << min_samples_overall << "\n"
        << "  Minor allele frequency cutoff: " << minor_allele_freq << "\n"
        << "  Maximum observed heterozygosity cutoff: " << max_obs_het << "\n"
        << "  Applying Fst correction: ";
    switch(fst_correction) {
    case p_value:
        fh << "P-value correction.\n";
        break;
    case bonferroni_win:
        fh << "Bonferroni correction within sliding window.\n";
        break;
    case bonferroni_gen:
        fh << "Bonferroni correction across genome wide sites.\n";
        break;
    case no_correction:
        fh << "none.\n";
        break;
    }
    fh << "  Pi/Fis kernel smoothing: " << (smooth_popstats == true ? "on" : "off") << "\n"
       << "  Fstats kernel smoothing: " << (smooth_fstats   == true ? "on" : "off") << "\n"
       << "  Bootstrap resampling: ";
    if (bootstrap)
        fh << "on, " << (bootstrap_type == bs_exact ? "exact; " : "approximate; ") << bootstrap_reps << " reptitions\n";
    else
        fh << "off\n";

    fh << "\n";
}

int

parse_command_line(int argc, char* argv[])
{
    bool   no_hap_exports = false;
    string map_type, map_format;
    char*  tmp_str;

    while (1) {
        static struct option long_options[] = {
            {"help",              no_argument,       NULL, 'h'},
            {"version",           no_argument,       NULL, 'v'},
            {"quiet",             no_argument,       NULL, 'q'},
            {"verbose",           no_argument,       NULL, 'd'},
            {"vcf",               no_argument,       NULL, 1004},
            {"fasta-loci",        no_argument,       NULL, 1006},
            {"fasta-samples",     no_argument,       NULL, 'J'},
            {"fasta-samples-raw", no_argument,       NULL, 'F'},
            {"structure",         no_argument,       NULL, 'S'},
            {"radpainter",        no_argument,       NULL, 1015}, {"fineRADstructure", no_argument, NULL, 1015},
            // {"fastphase",      no_argument,       NULL, 'A'},
            // {"phase",          no_argument,       NULL, 'C'},
            {"plink",          no_argument,       NULL, 'K'},
            {"genepop",        no_argument,       NULL, 1010},
            {"genepop-haps-3digits", no_argument, NULL, 1011},
            {"phylip",         no_argument,       NULL, 'Y'},
            {"phylip-var",     no_argument,       NULL, 'L'},
            {"phylip-var-all", no_argument,       NULL, 'T'},
            {"hzar",           no_argument,       NULL, 'Z'},
            {"treemix",        no_argument,       NULL, 'U'},
            {"merge-sites",    no_argument,       NULL, 'D'},
            {"sigma",          required_argument, NULL, 1005},
            {"threads",        required_argument, NULL, 't'},
            {"in-path",        required_argument, NULL, 'P'},
            {"v1",             no_argument,       NULL, 2000},
            {"out-path",       required_argument, NULL, 'O'},
            {"in-vcf",         required_argument, NULL, 'V'},
            {"min-samples-overall",   required_argument, NULL, 'R'},
            {"min-samples-per-pop",   required_argument, NULL, 'r'},
            {"min-populations",       required_argument, NULL, 'p'},
            {"filter-haplotype-wise", no_argument,       NULL, 'H'},
            {"renz",              required_argument, NULL, 'e'},
            {"popmap",            required_argument, NULL, 'M'},
            {"no-popmap",         no_argument,       NULL, 1017}, // Negates a previous -M/--popmap
            {"whitelist",         required_argument, NULL, 'W'},
            {"blacklist",         required_argument, NULL, 'B'},
            {"batch-size",        required_argument, NULL, 1999},
            {"dbg-min-gt-depth",  required_argument, NULL, 2001},
            {"write-single-snp",  no_argument,       NULL, 'I'}, 
            {"write-random-snp",  no_argument,       NULL, 'j'}, 
            {"no-hap-exports",    no_argument,       NULL, 1012}, 
            {"ordered-export",    no_argument,       NULL, 1002}, 
            {"smooth",            no_argument,       NULL, 'k'},
            {"smooth-fstats",     no_argument,       NULL, 1007}, 
            {"smooth-popstats",   no_argument,       NULL, 1008}, 
            {"fstats",            no_argument,       NULL, '6'},
            {"hwe",               no_argument,       NULL, 1014},
            {"log-fst-comp",      no_argument,       NULL, 'l'}, 
            {"bootstrap-type",    required_argument, NULL, 1001},
            {"bootstrap-reps",    required_argument, NULL, 1003},
            {"bootstrap-wl",      required_argument, NULL, 'Q'},
            {"bootstrap",         no_argument,       NULL, '1'},
            {"bootstrap-fst",     no_argument,       NULL, '2'},
            {"bootstrap-phist",   no_argument,       NULL, '3'},
            {"bootstrap-div",     no_argument,       NULL, '4'},
            {"bootstrap-pifis",   no_argument,       NULL, '5'},
            {"min-maf",           required_argument, NULL, 'a'},
            {"min-mac",           required_argument, NULL, 1016},
            {"max-obs-het",       required_argument, NULL, 1013},
            {"merge-prune-lim",   required_argument, NULL, 'i'},
            {"fst-correction",    no_argument,       NULL, 'f'},
            {"p-value-cutoff",    required_argument, NULL, 'u'},
            {"map-type",          required_argument, NULL, 3000},
            {"map-format",        required_argument, NULL, 3001},
            {"debug-flags",       required_argument, NULL, 1000},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int c = getopt_long(argc, argv, "ACDFHJKLNSTUV:YZ123456dhjklnqva:c:e:f:i:o:p:r:t:u:w:B:I:M:O:P:R:Q:W:", long_options, NULL);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'v':
            version();
            exit(1);
            break;
        case 'h':
            help();
            break;
        case 'q':
            quiet = true;
            break;
        case 'd':
            verbose = true;
            break;
        case 't':
            num_threads = atoi(optarg);
            break;
        case 'P':
            in_path = optarg;
            if (!in_path.empty() && in_path.back() != '/')
                in_path += "/";
            break;
        case 2000: //v1
            input_mode = InputMode::stacks;
            break;
        case 1999:
            batch_size = is_integer(optarg);
            break;
        case 2001:
            min_gt_depth = stol(optarg);
            break;
        case 'O':
            out_path = optarg;
            if (!out_path.empty() && out_path.back() != '/')
                out_path += "/";
            break;
        case 'V':
            in_vcf_path = optarg;
            break;
        case 'M': // popmap
            pmap_path = optarg;
            break;
        case 1017: //no-popmap
            pmap_path.clear();
            break;
        case 'D':
            merge_sites = true;
            break;
        case 'i':
            merge_prune_lim = is_double(optarg);
            if (merge_prune_lim > 1.0)
                merge_prune_lim = merge_prune_lim / 100;

            if (merge_prune_lim < 0 || merge_prune_lim > 1.0) {
                cerr << "Error: Unable to parse the merge sites pruning limit.\n";
                help();
            }
            break;
        case 1013:
            max_obs_het = is_double(optarg);
            if (max_obs_het > 1)
                max_obs_het = max_obs_het / 100;

            if (max_obs_het < 0 || max_obs_het > 1.0) {
                cerr << "Error: Unable to parse the maximum observed heterozygosity.\n";
                help();
            }
            break;
        case 'r':
            min_samples_per_pop = stod(optarg);
            if (!std::isfinite(min_samples_per_pop)) {
                cerr << "Error: Invalid value for --min-samples-per-pop, \""
                     << optarg << "\"\n";
                help();
            }
            if (min_samples_per_pop > 1.0)
                min_samples_per_pop /= 100.0;
            if (min_samples_per_pop < 0.0 || min_samples_per_pop > 1.0) {
                cerr << "Error: Invalid value for --min-samples-per-pop, \""
                     << optarg << "\"\n";
                help();
            }
            break;
        case 'R':
            min_samples_overall = stod(optarg);
            if (!std::isfinite(min_samples_overall)) {
                cerr << "Error: Invalid value for --min-samples-overall, \""
                     << optarg << "\"\n";
                help();
            }
            if (min_samples_overall > 1.0)
                min_samples_overall /= 100.0;
            if (min_samples_overall < 0.0 || min_samples_overall > 1.0) {
                cerr << "Error: Invalid value for --min-samples-overall, \""
                     << optarg << "\"\n";
                help();
            }
            break;
        case 'p':
            min_populations = stoi(optarg);
            break;
        case 'H':
            filter_haplotype_wise = true;
            break;
        case 'k':
            smooth_popstats = true;
            smooth_fstats   = true;
            calc_fstats     = true;
            ordered_export  = true;
            break;
        case 1007:
            smooth_fstats   = true;
            calc_fstats     = true;
            ordered_export  = true;
            break;
        case 1008:
            smooth_popstats = true;
            ordered_export  = true;
            break;
        case '6':
            calc_fstats = true;
            break;
        case 'l':
            log_fst_comp = true;
            break;
        case '1':
            ordered_export  = true;
            bootstrap       = true;
            bootstrap_fst   = true;
            bootstrap_phist = true;
            bootstrap_pifis = true;
            bootstrap_div   = true;
            break;
        case '2':
            ordered_export  = true;
            bootstrap_fst   = true;
            break;
        case '3':
            ordered_export  = true;
            bootstrap_phist = true;
            break;
        case '4':
            ordered_export  = true;
            bootstrap_div = true;
            break;
        case '5':
            ordered_export  = true;
            bootstrap_pifis = true;
            break;
        case 1001:
            if (strcasecmp(optarg, "exact") == 0)
                bootstrap_type = bs_exact;
            else if (strcasecmp(optarg, "approx") == 0)
                bootstrap_type = bs_approx;
            else {
                cerr << "Error: Unknown bootstrap type specified '" << optarg << "'\n";
                help();
            }
            break;
        case 1003:
            ordered_export  = true;
            bootstrap = true;
            bootstrap_reps = atoi(optarg);
            break;
        case 'Q':
            bs_wl_file = optarg;
            bootstrap_wl = true;
            break;
        case 'I':
            write_single_snp = true;
            no_hap_exports = true;
            break;
        case 'j':
            write_random_snp = true;
            no_hap_exports = true;
            break;
        case 1012: // --no-haps
            no_hap_exports = true;
            break;
        case 1002:
            ordered_export = true;
            break;
        case 1004: // --vcf
            add_export<VcfExport>();
            add_export<VcfHapsExport>();
            break;
        case 1006: // --fasta-loci
            add_export<FastaLociExport>();
            break;
        case 'F':
            add_export<FastaRawExport>();
            break;
        case 'J':
            add_export<FastaSamplesExport>();
            break;
        case 1010: // --genepop
            add_export<GenePopExport>();
            add_export<GenePopHapsExport>();
            break;
        case 1011: //genepop-haps-3digits
            add_export<GenePopHapsExport>();
            dynamic_cast<GenePopHapsExport&>(**find_export<GenePopHapsExport>()).set_digits(3);
            break;
        case 1015: // --radpainter,--fineradstructure
            add_export<FineRADStructureExport>();
            break;
        case 1014: // --hwe
            calc_hwp = true;
            break;
        case 'S':
            add_export<StructureExport>();
            break;
        case 'A':
            cerr << "Ignoring --fastphase output request, which is not currently implemented.\n";
            break;
        case 'C':
            cerr << "Ignoring --phase output request, which is not currently implemented.\n";
            break;
        case 'K': // --plink
            add_export<PlinkExport>();
            break;
        case 'Z':
            add_export<HzarExport>();
            break;
        case 'Y':
            add_export<PhylipFixedExport>();
            break;
        case 'L':
            add_export<PhylipVarExport>();
            break;
        case 'T':
            add_export<PhylipVarAllExport>();
            break;
        case 'U':
            add_export<TreemixExport>();
            break;
        case 'W':
            wl_file = optarg;
            break;
        case 'B':
            bl_file = optarg;
            break;
        case 3000:
            map_type = optarg;
            break;
        case 3001:
            map_format = optarg;
            break;
        case 'a':
            minor_allele_freq = atof(optarg);
            if (minor_allele_freq > 1)
                minor_allele_freq = minor_allele_freq / 100;

            if (minor_allele_freq < 0 || minor_allele_freq > 0.5) {
                cerr << "Error: Unable to parse the minor allele frequency.\n";
                help();
            }
            break;
        case 1016:
            minor_allele_cnt = strtol(optarg, &tmp_str, 10);
            if (tmp_str == optarg || *tmp_str != '\0' || minor_allele_cnt < 0) {
                cerr << "Error: Unable to parse the minor allele count from '" << optarg << "'.\n";
                help();
            }
            break;
        case 'f':
            fst_correction = p_value;
            // if (strcasecmp(optarg, "p_value") == 0)
            //     fst_correction = p_value;
            // else if (strcasecmp(optarg, "bonferroni_win") == 0)
            //     fst_correction = bonferroni_win;
            // else if (strcasecmp(optarg, "bonferroni_gen") == 0)
            //     fst_correction = bonferroni_gen;
            // else {
            //     cerr << "Error: Unknown Fst correction specified '" << optarg << "'\n";
            //     help();
            // }
            break;
        case 'u':
            p_value_cutoff = atof(optarg);
            break;
        case 'e':
            enz = optarg;
            enz.at(0) = tolower(enz.at(0));
            if (renz.count(enz) == 0) {
                cerr << "Error: Unrecognized restriction enzyme specified: '" << enz.c_str() << "'.\n";
                help();
            }
            break;
        case 1005: //sigma
            sigma = atof(optarg);
            break;
        case '?':
            // getopt_long already printed an error message.
            help();
            break;
        case 1000:
        {
            static const set<string> known_debug_flags = {"VCFCOMP"};
            stringstream ss (optarg);
            string s;
            while (getline(ss, s, ',')) {
                if (known_debug_flags.count(s)) {
                    debug_flags.insert(s);
                } else {
                    cerr << "DEBUG: Error: Unknown error flag '" << s << "'.\n";
                    return -1;
                }
            }
            cout << "DEBUG: Debug flag(s) : '" << optarg << "'.\n";

            if (debug_flags.count("VCFCOMP") && not write_random_snp) {
                write_single_snp = true;
                cout << "DEBUG: Added --write-single-snp.\n";
            }

            break;
        }
        default:
            cerr << "Error: Unknown command line option: '" << (char) c << "'\n";
            help();
            exit(1);
        }
    }

    if (optind < argc) {
        cerr << "Error: Failed to parse command line: '" << argv[optind] << "' is seen as a positional argument. Expected no positional arguments.\n";
        help();
    }

    //
    // Check argument constraints.
    //
    if (!in_path.empty() && !in_vcf_path.empty()) {
        cerr << "Error: Please specify either '-P/--in-path' or '-V/--in-vcf', not both.\n";
        help();
    } else if (in_path.empty() && in_vcf_path.empty()) {
        cerr << "Error: One of '-P/--in-path' or '-V/--in-vcf' is required.\n";
        help();
    } else if (not in_vcf_path.empty()) {
        input_mode = InputMode::vcf;
    }

    if (input_mode == InputMode::stacks || input_mode == InputMode::stacks2) {
        if (out_path.empty())
            out_path = in_path;

        out_prefix = "populations";

    } else if (input_mode == InputMode::vcf) {

        if (out_path.empty()) {
            cerr << "Error: Malformed arguments: input mode 'vcf' requires an output directory (--out-path).\n";
            help();
        }

        // Determine out_prefix
        string fname = in_vcf_path;
        if (in_vcf_path.find_last_of('/') != string::npos && in_vcf_path.back() != '/')
            fname = in_vcf_path.substr(in_vcf_path.find_last_of('/')+1);
        size_t trim = 0;
        if (fname.length() > 4 && fname.substr(fname.length()-4) == ".vcf")
            trim = 4;
        else if (fname.length() > 7 && fname.substr(fname.length()-7) == ".vcf.gz")
            trim = 7;
        out_prefix = fname.substr(0, fname.length()-trim);
        out_prefix += ".p";
    }

    // Other
    if (write_single_snp + write_random_snp + filter_haplotype_wise > 1) {
        cerr << "Error: At most one of --write-single-snp, --write-random-snp"
             << " or --filter-haplotype-wise may be specified.\n";
        help();
    }

    //
    // Genetic map export option constraints.
    //
    if ((map_type.empty() && !map_format.empty()) || (!map_type.empty() && map_format.empty())) {
        cerr << "Error: for genetic map export you must specify both the map type (--map-type) and format (--map-format).\n";
        help();
    }
    if (!map_type.empty() && !map_format.empty()) {
        mapping_cross = true;

        transform(map_type.begin(), map_type.end(), map_type.begin(), ::tolower);

        if (map_type == "cp")
            mapcross_type = CrossT::cp;
        else if (map_type == "bc1")
            mapcross_type = CrossT::bc1;
        else if (map_type == "f2")
            mapcross_type = CrossT::f2;
        else if (map_type == "dh")
            mapcross_type = CrossT::dh;
        else {
            cerr << "Error: Unrecognized mapping cross type specified, '" << map_type << "'.\n";
            help();
        }

        transform(map_format.begin(), map_format.end(), map_format.begin(), ::tolower);

        if (map_format == "joinmap")
            mapcross_format = FormatT::joinmap;
        else if (map_format == "rqtl")
            mapcross_format = FormatT::rqtl;
        else if (map_format == "onemap")
            mapcross_format = FormatT::onemap;
        else {
            cerr << "Error: Unrecognized mapping cross output format specified, '" << map_format << "'.\n";
            help();
        }
    }

    if (merge_sites == true && enz.length() == 0) {
        cerr << "Error: You must specify the restriction enzyme associated with this data set to merge overlaping cutsites.\n";
        help();
    }

    if (no_hap_exports) {
        for (Export*& e : exports) {
            if (e->is_hap_export()) {
                delete e;
                e = NULL;
            }
        }
        stacks_erase_if(exports, [](const Export* e){return e == NULL;});
    }

    return 0;
}

void version() {
    cout << "populations " << VERSION << "\n\n";
}

void help() {
    cout << "populations " << VERSION << "\n"
         << "Usage:\n"
         << "populations -P dir [-O dir] [-M popmap] (filters) [--fstats] [-k [--sigma=150000] [--bootstrap [-N 100]]] (output formats)\n"
         << "populations -V vcf -O dir [-M popmap] (filters) [--fstats] [-k [--sigma=150000] [--bootstrap [-N 100]]] (output formats)\n"
         << "\n"
         << "  -P,--in-path: path to a directory containing Stacks ouput files.\n"
         << "  -V,--in-vcf: path to a standalone input VCF file.\n"
         << "  -O,--out-path: path to a directory where to write the output files. (Required by -V; otherwise defaults to value of -P.)\n"
         << "  -M,--popmap: path to a population map. (Format is 'SAMPLE1 \\t POP1 \\n SAMPLE2 ...'.)\n"
         << "  -t,--threads: number of threads to run in parallel sections of code.\n"
         << "  --batch-size [int]: the number of loci to process in a batch (default: 10,000 in de novo mode; in reference mode, one chromosome\n"
         << "                      per batch). Increase to speed analysis, uses more memory, decrease to save memory).\n"
         << "\n"
         << "Data Filtering:\n"
         << "  -p,--min-populations [int]: minimum number of populations a locus must be present in to process a locus.\n"
         << "  -r,--min-samples-per-pop [float]: minimum percentage of individuals in a population required to process a locus for that population.\n"
         << "  -R,--min-samples-overall [float]: minimum percentage of individuals across populations required to process a locus.\n"
         << "  -H,--filter-haplotype-wise: apply the above filters haplotype wise\n"
         << "                              (unshared SNPs will be pruned to reduce haplotype-wise missing data).\n"
         << "  --min-maf [float]: specify a minimum minor allele frequency required to process a SNP (0 < min_maf < 0.5).\n"
         << "  --min-mac [int]: specify a minimum minor allele count required to process a SNP.\n"
         << "  --max-obs-het [float]: specify a maximum observed heterozygosity required to process a SNP.\n"
         << "\n"
         << "  --write-single-snp: restrict data analysis to only the first SNP per locus (implies --no-haps).\n"
         << "  --write-random-snp: restrict data analysis to one random SNP per locus (implies --no-haps).\n"
         << "  -B,--blacklist: path to a file containing Blacklisted markers to be excluded from the export.\n"
         << "  -W,--whitelist: path to a file containing Whitelisted markers to include in the export.\n"
         << "\n"
         << "Merging and Phasing:\n"
         << "  -e,--renz: restriction enzyme name.\n"
         << "  --merge-sites: merge loci that were produced from the same restriction enzyme cutsite (requires reference-aligned data).\n"
         << "  --merge-prune-lim: when merging adjacent loci, if at least X% samples posses both loci prune the remaining samples out of the analysis.\n"
         << "\n"
         << "Locus stats:\n"
         << "  --hwe: calculate divergence from Hardy-Weinberg equilibrium for each locus.\n"
         << "\n"
         << "Fstats:\n"
         << "  --fstats: enable SNP and haplotype-based F statistics.\n"
         << "  --fst-correction: specify a p-value correction to be applied to Fst values based on a Fisher's exact test. Default: off.\n"
         << "  --p-value-cutoff [float]: maximum p-value to keep an Fst measurement. Default: 0.05. (Also used as base for Bonferroni correction.)\n"
         << "\n"
         << "Kernel-smoothing algorithm:\n"
         << "  -k,--smooth: enable kernel-smoothed Pi, Fis, Fst, Fst', and Phi_st calculations.\n"
         << "  --smooth-fstats: enable kernel-smoothed Fst, Fst', and Phi_st calculations.\n"
         << "  --smooth-popstats: enable kernel-smoothed Pi and Fis calculations.\n"
         << "    (Note: turning on smoothing implies --ordered-export.)\n"
         << "  --sigma [int]: standard deviation of the kernel smoothing weight distribution. Default 150kb.\n"
         << "  --bootstrap: turn on boostrap resampling for all smoothed statistics.\n"
         << "  -N,--bootstrap-reps [int]: number of bootstrap resamplings to calculate (default 100).\n"
         << "  --bootstrap-pifis: turn on boostrap resampling for smoothed SNP-based Pi and Fis calculations.\n"
         << "  --bootstrap-fst: turn on boostrap resampling for smoothed Fst calculations based on pairwise population comparison of SNPs.\n"
         << "  --bootstrap-div: turn on boostrap resampling for smoothed haplotype diveristy and gene diversity calculations based on haplotypes.\n"
         << "  --bootstrap-phist: turn on boostrap resampling for smoothed Phi_st calculations based on haplotypes.\n"
         << "  --bootstrap-wl [path]: only bootstrap loci contained in this whitelist.\n"
         << "\n"
         << "File output options:\n"
         << "  --ordered-export: if data is reference aligned, exports will be ordered; only a single representative of each overlapping site.\n"
         << "  --fasta-loci: output locus consensus sequences in FASTA format.\n"
         << "  --fasta-samples: output the sequences of the two haplotypes of each (diploid) sample, for each locus, in FASTA format.\n"
         << "  --vcf: output SNPs and haplotypes in Variant Call Format (VCF).\n"
         << "  --genepop: output SNPs and haplotypes in GenePop format.\n"
         << "  --structure: output results in Structure format.\n"
         << "  --radpainter: output results in fineRADstructure/RADpainter format.\n"
         // << "  --phase*: output genotypes in PHASE format.\n"
         // << "  --fastphase*: output genotypes in fastPHASE format.\n"
         << "  --plink: output genotypes in PLINK format.\n"
         << "  --hzar: output genotypes in Hybrid Zone Analysis using R (HZAR) format.\n"
         << "  --phylip: output nucleotides that are fixed-within, and variant among populations in Phylip format for phylogenetic tree construction.\n"
         << "  --phylip-var: include variable sites in the phylip output encoded using IUPAC notation.\n"
         << "  --phylip-var-all: include all sequence as well as variable sites in the phylip output encoded using IUPAC notation.\n"
         << "  --treemix: output SNPs in a format useable for the TreeMix program (Pickrell and Pritchard).\n"
         << "  --no-hap-exports: omit haplotype outputs.\n"
         << "  --fasta-samples-raw: output all haplotypes observed in each sample, for each locus, in FASTA format.\n"
         << "\n"
         << "Genetic map output options (population map must specify a genetic cross):\n"
         << "  --map-type: genetic map type to write. 'CP', 'DH', 'F2', and 'BC1' are the currently supported map types.\n"
         << "  --map-format: mapping program format to write, 'joinmap', 'onemap', and 'rqtl' are currently supported.\n"
         << "\n"
         << "Additional options:\n"
         << "  -h,--help: display this help messsage.\n"
         << "  -v,--version: print program version.\n"
         << "  --verbose: turn on additional logging.\n"
         << ("  --log-fst-comp: log components of Fst/Phi_st calculations to a file.\n")
         #ifdef DEBUG
         << "\n"
         << "DEBUG:\n"
         << "  --dbg-min-gt-depth\n"
         << "  --genepop-haps-3digits: Use 3-digit alleles in the genepop haps output (default is 2-digit).\n"
         #endif
         ;

    // << "    --bootstrap-type [exact|approx]: enable bootstrap resampling for population statistics (reference genome required).\n"

    exit(1);
}
