// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2019-2020, Julian Catchen <jcatchen@illinois.edu>
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
// mapping_utils.cc -- common routines for genetic mapping crosses.
//

#include "mapping_utils.h"

MappingGenotypesProcessor::MappingGenotypesProcessor(MetaPopInfo *popi, CrossT ct)
{
    this->_mpopi      = popi;
    this->_cross_type = ct;

    this->_num_loci         = 0;
    this->_mappable_loci    = 0;
    this->_mean_progeny_sum = 0.0;

    switch(this->_cross_type) {
    case CrossT::dh:
        load_dh_dictionary(this->_marker_map, this->_genotype_map);
        break;
    case CrossT::cp:
        load_cp_dictionary(this->_marker_map, this->_genotype_map);
        break;
    case CrossT::bc1:
        load_bc_dictionary(this->_marker_map, this->_genotype_map);
        break;
    case CrossT::f2:
        load_f2_dictionary(this->_marker_map, this->_genotype_map);
        break;
    case CrossT::unk:
    default:
        break;
    }

    //
    // Load segregation ratios for this map type (for calculating segregation distortion).
    //
    load_segregation_ratios(this->_cross_type, this->_segregation_ratios);
    
    //
    // Find the 'progeny' population.
    //
    const vector<Pop> &pops = this->_mpopi->pops();
    uint index;
    for (index = 0; index < pops.size(); index++)
        if (pops[index].name == "progeny")
            break;
    this->_progeny_index = index;
}

bool
MappingGenotypesProcessor::is_mapping_cross(ostream &log_fh)
{
    //
    // Assess whether the population map only specifies parents (no more than 2) and progeny.
    //
    const vector<Pop> &pops = this->_mpopi->pops();
    if (pops.size() != 2 || (pops[0].name != "parent" && pops[1].name != "parent")) {
        log_fh << "Error: for genetic map export population map must specify only two groupings, 'parent' and 'progeny'.\n";
        return false;
    }

    bool progeny_found = false;
    
    for (uint i = 0; i < pops.size(); i++) {
        if (pops[i].name == "parent" && pops[i].n_samples() > 2) {
            log_fh << "Error: for genetic map export there can be no more than two 'parent' samples.\n";
            return false;
        } else if (pops[i].name == "progeny")
            progeny_found = true;
    }

    return progeny_found;
}

int
MappingGenotypesProcessor::next_batch(const vector<LocBin *> &loci)
{
    size_t fsample = this->_mpopi->pops()[this->_progeny_index].first_sample;
    size_t lsample = this->_mpopi->pops()[this->_progeny_index].last_sample;

    for (uint k = 0; k < loci.size(); k++) {
        const LocBin  *loc  = loci[k];
        CSLocus      *cloc = loc->cloc;

        this->_num_loci++;

        //
        // Classify this marker.
        //
        string marker = classify_marker(loc, this->_mpopi);

        if (marker.length() == 0 || this->_marker_map.count(marker) == 0)
            continue;

        this->_mappable_loci++;

        map<string, string> parent_gmap;
        assign_parent_genotypes(this->_mpopi, loc, marker, parent_gmap);

        cloc->marker = this->_marker_map[marker];

        Datum **d = loc->d;

        map<string, int> genotype_cnts;

        //
        // Assign genotypes to the mapping population.
        //
        for (uint j = fsample; j <= lsample; j++) {
            if (d[j] == NULL)
                continue;
            
            string g = assign_generic_genotype(marker, d[j], parent_gmap);
            d[j]->genotype(g);
            cloc->cnt++;
            
            if (g[0] != '-') {
                genotype_cnts[g]++;
                cloc->gcnt++;
                this->_mean_progeny_sum++;
            }
        }

        //
        // For CP maps, use segregation ratios to detect and correct cases of missing parental alleels. 
        //
        if (this->_cross_type == CrossT::cp)
            this->_corrected += this->correct_cp_marker(cloc, d, fsample, lsample, genotype_cnts);
        
        //
        // Quantify the amount of segregation distortion from Hardy-Weinberg equilibrium using a chi-squared test.
        //
        cloc->chisq = chisq_test(this->_segregation_ratios, genotype_cnts, cloc->marker, cloc->gcnt);
    }

    return 0;
}

string
classify_marker(const LocBin *loc, const MetaPopInfo *mpopi)
{
    const Datum *d_1, *d_2;

    //
    // Count the number of parental tags matching this catalog tag. A proper marker should
    // contain a single representative from each parent; multiple alleles must be called from
    // a single tag from a single parent.
    //
    const vector<Pop> &pops = mpopi->pops();
    uint index;
    for (index = 0; index < pops.size(); index++)
        if (pops[index].name == "parent")
            break;

    size_t parent_count = pops[index].n_samples();
    if (parent_count == 1) {
        d_1 = loc->d[pops[index].first_sample];
        d_2 = NULL;
    } else {
        d_1 = loc->d[pops[index].first_sample];
        d_2 = loc->d[pops[index].last_sample];
    }

    //
    // How many parents have a genotype for this locus?
    //
    parent_count = 0;
    if (d_1 != NULL) parent_count++;
    if (d_2 != NULL) parent_count++;

    string marker;

    if (parent_count == 2) {
        //
        // Locus is present in both parents.
        //

        //
        // Determine the number of unique alleles in each parent and overall.
        //
        set<string> unique_alleles_1, unique_alleles_2, unique_alleles;

        for (auto hit = d_1->obshap.begin(); hit != d_1->obshap.end(); hit++)
            unique_alleles_1.insert(*hit);
        for (auto hit = d_2->obshap.begin(); hit != d_2->obshap.end(); hit++)
            unique_alleles_2.insert(*hit);
        unique_alleles.insert(unique_alleles_1.begin(), unique_alleles_1.end());
        unique_alleles.insert(unique_alleles_2.begin(), unique_alleles_2.end());
        size_t allele_cnt_1       = unique_alleles_1.size();
        size_t allele_cnt_2       = unique_alleles_2.size();
        size_t num_unique_alleles = unique_alleles.size();

        if (allele_cnt_1 == 2 && allele_cnt_2 == 2) {
            //
            // Locus is heterozygous in both parents. However, the number of alleles present distinguishes
            // what type of marker it is. Four unique alleles requries an ab/cd marker, while four
            // alleles that are the same in both parents requires an ab/ab marker. Finally, three unique
            // alleles requires either an ab/ac marker.
            //
            if (num_unique_alleles == 3)
                marker = "ab/ac";
            else if (num_unique_alleles == 2)
                marker = "ab/ab";
            else
                marker = "ab/cd";

        } else if (allele_cnt_1 == 2 && allele_cnt_2 == 1) {
            //
            // Locus is homozygous in one parent and heterozygous in the other.
            //
            if (num_unique_alleles == 3)
                marker = "ab/cc";
            else if (num_unique_alleles == 2)
                marker = "ab/aa";

        } else if (allele_cnt_1 == 1 && allele_cnt_2 == 2) {
            //
            // Locus is homozygous in one parent and heterozygous in the other.
            //
            if (num_unique_alleles == 3)
                marker = "cc/ab";
            else if (num_unique_alleles == 2)
                marker = "aa/ab";

        } else if (allele_cnt_1 == 1 && allele_cnt_2 == 1) {
            //
            // Locus is homozygous in both parents, but heterozygous between parents.
            //
            if (strcmp(d_1->obshap[0], d_2->obshap[0]) != 0)
                marker = "aa/bb";
        }

    } else if (parent_count == 1) {
        //
        // Locus only exists in one parent.
        //
        if (d_1 != NULL && d_1->obshap.size() == 2)
            marker = "ab/--";
        else if (d_2 != NULL && d_2->obshap.size() == 2)
            marker = "--/ab";
    }

    return marker;
}

int
assign_parent_genotypes(const MetaPopInfo *mpopi,
                        const LocBin *loc,
                        string marker,
                        map<string, string> &gmap)
{
    //
    // Create a genotype map. For any set of alleles, this routine will
    // assign each allele to one of the constituent genotypes, e.g. given the
    // marker type 'aaxbb' and the alleles 'A' from the male, and 'G'
    // from the female, will assign 'G' == 'bb' and 'A'== 'aa'. It assumes that
    // recombination may have occurred as with an F2, F3 or later cross.
    //

    //
    // Identify the parents.
    //
    Datum *d_1, *d_2;
    uint   index;
    const vector<Pop> &pops = mpopi->pops();
    for (index = 0; index < pops.size(); index++)
        if (pops[index].name == "parent")
            break;

    d_1 = loc->d[pops[index].first_sample];
    d_2 = loc->d[pops[index].last_sample];
    
    set<char> p1_gtypes, p2_gtypes;
    map<char, int> legal_gtypes, com_gtypes;

    assert(marker.length() == 5);

    //
    // First, identify any alleles that are common between the two parents.
    //
    p1_gtypes.insert(marker[0]);
    p1_gtypes.insert(marker[1]);
    p2_gtypes.insert(marker[3]);
    p2_gtypes.insert(marker[4]);

    for (auto i = p1_gtypes.begin(); i != p1_gtypes.end(); i++)
        if (*i != '-') legal_gtypes[*i]++;
    for (auto i = p2_gtypes.begin(); i != p2_gtypes.end(); i++)
        if (*i != '-') legal_gtypes[*i]++;
    //
    // Find the common genotypes
    //
    vector<char> types;
    for (auto j = legal_gtypes.begin(); j != legal_gtypes.end(); j++)
        if (j->second > 1) types.push_back(j->first);
    sort(types.begin(), types.end());

    set<string> p1_obshap, p2_obshap;
    map<string, int> haplotypes;
    vector<pair<string, int> > sorted_haplotypes;

    if (d_1 != NULL) {
        for (uint n = 0; n < d_1->obshap.size(); n++)
            p1_obshap.insert(d_1->obshap[n]);
        for (auto n = p1_obshap.begin(); n != p1_obshap.end(); n++)
            haplotypes[*n]++;
    }
    if (d_2 != NULL) {
        for (uint n = 0; n < d_2->obshap.size(); n++)
            p2_obshap.insert(d_2->obshap[n]);
        for (auto n = p2_obshap.begin(); n != p2_obshap.end(); n++)
            haplotypes[*n]++;
    }

    //
    // Sort the haplotypes map by value
    //
    for (auto k = haplotypes.begin(); k != haplotypes.end(); k++)
        sorted_haplotypes.push_back(*k);
    sort(sorted_haplotypes.begin(), sorted_haplotypes.end(),
         [](pair<string, int> a, pair<string, int> b) {
             return (a.second > b.second);
         });

    //
    // Identify any common alleles.
    //
    for (uint n = 0, index = 0; n < sorted_haplotypes.size() && index < types.size(); n++, index++) {
        if (sorted_haplotypes[n].second > 1) {
            gmap[sorted_haplotypes[n].first] = types[index];
            com_gtypes[types[index]]++;
            // cerr << "  Assigning common allele " << sorted_haplotypes[n].first << " to genotype '" << gmap[sorted_haplotypes[n].first] << "'\n";
        }
    }

    //
    // Now, examine the remaining first parent alleles.
    //
    if (d_1 != NULL) {
        legal_gtypes.clear();
        for (auto i = p1_gtypes.begin(); i != p1_gtypes.end(); i++)
            if (*i != '-' && com_gtypes.count(*i) == 0) {
                // cerr << "  Adding " << *i << " to first parent genotypes\n";
                legal_gtypes[*i]++;
            }
        types.clear();
        for (auto j = legal_gtypes.begin(); j != legal_gtypes.end(); j++)
            types.push_back(j->first);
        sort(types.begin(), types.end());

        for (uint n = 0, index = 0; n < d_1->obshap.size() && index < types.size(); n++, index++) {
            if (gmap.count(d_1->obshap[n])) {
                index--;
                continue;
            }
            gmap[d_1->obshap[n]] = types[index];
            // cerr << "  Assinging '" << d_1->obshap[n] << "' to first parent genotype '" << gmap[d_1->obshap[n]] << "'\n";
        }
    }

    //
    // Finally, repeat in the second parent.
    //
    if (d_2 != NULL) {
        legal_gtypes.clear();
        for (auto i = p2_gtypes.begin(); i != p2_gtypes.end(); i++)
            if (*i != '-' && com_gtypes.count(*i) == 0) {
                // cerr << "  Adding " << *i << " to second genotypes\n";
                legal_gtypes[*i]++;
            }
        types.clear();
        for (auto j = legal_gtypes.begin(); j != legal_gtypes.end(); j++)
            types.push_back(j->first);
        sort(types.begin(), types.end());

        for (uint n = 0, index = 0; n < d_2->obshap.size() && index < types.size(); n++, index++) {
            if (gmap.count(d_2->obshap[n])) {
                index--;
                continue;
            }
            gmap[d_2->obshap[n]] = types[index];
            // cerr << "  Assinging '" << d_2->obshap[n] << "' to second parent genotype '" << gmap[d_2->obshap[n]] << "'\n";
        }
    }

    return 0;
}

string
assign_generic_genotype(string marker,
                        const Datum *d,
                        map<string, string> &parent_gmap)
{
    string m;

    if (d == NULL)
        return m;

    vector<string> gtypes;
    string gtype;

    for (uint j = 0; j < d->obshap.size(); j++) {
        //
        // Impossible allele encountered.
        //
        if (parent_gmap.count(d->obshap[j]) == 0) {
            gtypes.clear();
            gtypes.push_back("-");
            break;
        }

        gtypes.push_back(parent_gmap[d->obshap[j]]);
    }

    sort(gtypes.begin(), gtypes.end());
    for (uint j = 0; j < gtypes.size(); j++)
        gtype += gtypes[j];

    return gtype;
}

double
chisq_test(map<string, map<string, double>> &seg_ratios, map<string, int> &cnts, string marker, double n)
{
    if (seg_ratios.count(marker) == 0)
        return 1.0;

    //
    // Calculate chi-square value.
    //   sit->second * n == the expected value for this genotype
    //
    double chisq  = 0.0;
    double exp    = 0.0;
    double obs    = 0.0;
    double df     = seg_ratios[marker].size() - 1;

    map<string, double>::iterator sit;

    for (sit = seg_ratios[marker].begin(); sit != seg_ratios[marker].end(); sit++) {
        obs = cnts.count(sit->first) == 0 ? 0 : cnts[sit->first];
        exp = sit->second * n;
        // cerr << "      category: " << sit->first << "; obs: " << obs << "; exp: " << exp << "\n";

        chisq += ((obs - exp) * (obs - exp)) / exp;
    }
    // cerr << "    df: " << df << "; Chisq value: " << chisq << "; pvalue: " << chisq_pvalue(df, chisq) << "\n";

    //
    // Determine p-value
    //
    return chisq_pvalue(df, chisq);
}

double
chisq_pvalue(int df, double chisq)
{
    int i = 0;
    while (chisq > chisq_crit_values[df][i] &&
           i < chisq_crit_values_size) {
        i++;
    }

    if (i == chisq_crit_values_size)
        return chisq_crit_values[0][chisq_crit_values_size - 1];

    return chisq_crit_values[0][i];
}

int
MappingGenotypesProcessor::correct_cp_marker(CSLocus *cloc, Datum **d, size_t fsample, size_t lsample, map<string, int> &genotype_cnts)
{
    if (cloc->gcnt == 0)
        return 0;
    
    //
    // When one of the parents of the cross is missing an allele, we
    // will see the ratios defined in this->_c_ratios_1. 
    //
    if (cloc->marker == "ab/aa" ||
        cloc->marker == "aa/ab" ||
        cloc->marker == "ab/cc" ||
        cloc->marker == "cc/ab") {

        //
        // Calculate initial segregation distortion.
        //
        double chisq_pval = chisq_test(this->_c_ratios_1, genotype_cnts, cloc->marker, cloc->gcnt);

        //
        // Check if our genotype ratios match the segregation ratios specified above. If so,
        // we have a dropped allele in one of the parents.
        //
        if (chisq_pval < chisq_pval_limit) {

            if (cloc->marker == "ab/aa")
                cloc->marker = "ab/a-";
            else if (cloc->marker == "aa/ab")
                cloc->marker = "-a/ab";
            else if (cloc->marker == "ab/cc")
                cloc->marker = "ab/c-";
            else if (cloc->marker == "cc/ab")
                cloc->marker = "-c/ab";

            if (cloc->marker == "ab/a-" || cloc->marker == "-a/ab") {

                for (uint i = fsample; i <= lsample; i++) {
                    if (d[i] == NULL) continue;

                    if (strcmp(d[i]->gtype, "bb") == 0)
                        strcpy(d[i]->gtype, "ab");
                }

            } else if (cloc->marker == "ab/c-") {

                for (uint i = fsample; i <= lsample; i++) {
                    if (d[i] == NULL) continue;

                    if (strcmp(d[i]->gtype, "bb") == 0)
                        strcpy(d[i]->gtype, "bd");
                    else if (strcmp(d[i]->gtype, "aa") == 0)
                        strcpy(d[i]->gtype, "ad");
                }

            } else if (cloc->marker == "-c/ab") {

                for (uint i = fsample; i <= lsample; i++) {
                    if (d[i] == NULL) continue;

                    if (strcmp(d[i]->gtype, "bb") == 0)
                        strcpy(d[i]->gtype, "ad");
                    else if (strcmp(d[i]->gtype, "aa") == 0)
                        strcpy(d[i]->gtype, "ac");
                    else if (strcmp(d[i]->gtype, "bc") == 0)
                        strcpy(d[i]->gtype, "bd");
                    else if (strcmp(d[i]->gtype, "ac") == 0)
                        strcpy(d[i]->gtype, "bc");
                }
            }
        }

        return 1;

    } else if (cloc->marker =="aa/bb") {
        //
        // Now we will deal with aa/bb markers separately, since there can be three possible
        // missing allele situations:
        //   aa: 50%, ab: 50% - we have an aa/b- marker, which should be mapped as an --/ab
        //   bb: 50%, ab: 50% - we have an -a/bb marker, which should be mapped as an ab/--
        //   aa: 33%, ab: 33%, bb: 33% - we have an -a/b- maker, which should be mapped as an ab/ab, but
        //                               we can't disambiguate the aa bb genotypes so it can't be mapped.
        //
        double chisq_pval = chisq_test(this->_c_ratios_2, genotype_cnts, cloc->marker, cloc->gcnt);

        if (chisq_pval >= chisq_pval_limit) {

            cloc->marker = "aa/b-";

            for (uint i = fsample; i <= lsample; i++) {
                if (d[i] == NULL) continue;

                if (strcmp(d[i]->gtype, "ab") == 0)
                    strcpy(d[i]->gtype, "bb");
            }
            
            return 1;

        } else {
            chisq_pval = chisq_test(this->_c_ratios_3, genotype_cnts, cloc->marker, cloc->gcnt);

            if (chisq_pval >= chisq_pval_limit) {
                cloc->marker = "-a/bb";

                for (uint i = fsample; i <= lsample; i++) {
                    if (d[i] == NULL) continue;

                    if (strcmp(d[i]->gtype, "ab") == 0)
                        strcpy(d[i]->gtype, "aa");
                }

                return 1;
            }
        }
    }

    return 0;
}
