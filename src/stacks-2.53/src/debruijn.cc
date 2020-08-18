#include <iostream>
#include <iomanip>
#include <sstream>

#include "constants.h"
#include "debruijn.h"

using namespace std;

SPath::SPath(Node* first) : first_(first), last_(NULL), d_(), visitdata(NULL) {
    d_.n_nodes = 1;
    d_.km_cumcount = first->count();
    Node* n = first_;
    Node* s = n->first_succ();
    while (n->n_succ() == 1 && s->n_pred() == 1 && s != first_) {
        // Extend the simple path.
        n = s;
        ++d_.n_nodes;
        d_.km_cumcount += n->count();
        s = n->first_succ();
    }
    last_ = n;
}

void SPath::erase(size_t km_len) {
    assert(!empty());
    if (first_ != last_) {
        assert(first_->n_succ() == 1);
        assert(last_->n_pred() == 1);
    }
    // Disconnect the SPath from the graph. Record neighbors.
    array<SPath*,4> preds{}, succs{};
    for (size_t nt2=0; nt2<4; ++nt2) {
        SPath* p = pred(nt2);
        SPath* s = succ(nt2);
        if (p != NULL && p != this) {
            preds[nt2] = p;
            p->last_->rm_succ(size_t(first_->km().back(km_len)), first_);
        }
        if (s != NULL && s != this) {
            succs[nt2] = s;
            last_->rm_succ(nt2, s->first_);
        }
    }
    // Clear the SPath.
    Node* next = NULL;
    for (Node* n = first_; ; n = next) {
        if (n != last_) {
            assert(n->n_succ() == 1);
            next = n->first_succ();
        }
        *n = Node();
        if (n == last_)
            break;
    }
    *this = SPath();
    // Finally, check if we need to merge some SPaths (i.e. if we created some
    // new "simple edges").
    for (size_t nt2=0; nt2<4; ++nt2) {
        SPath* p = preds[nt2];
        if (p != NULL && !p->empty())
            p->merge_forward();
    }
    for (size_t nt2=0; nt2<4; ++nt2) {
        SPath* s = succs[nt2];
        if (s != NULL && !s->empty() && s->n_pred() == 1) {
            SPath* s_pred = NULL;
            for (size_t nt22=0; nt22<4; ++nt22) {
                if (s->pred(nt22) != NULL) {
                    s_pred = s->pred(nt22);
                    break;
                }
            }
            assert(s_pred != NULL && !s_pred->empty());
            s_pred->merge_forward();
        }
    }
}

bool SPath::merge_forward() {
    assert(!empty());
    if (!is_mergeable_forward())
        return false;
    // Clear the SPaths we can merge into. (Keep track of the predicted end.)
    SPath* s = first_succ();
    assert(s != this);
    Node* end_node = s->last_;
    while (s->is_mergeable_forward()) {
        SPath* tmp = s->first_succ();
        if (tmp == this)
            break;
        *s = SPath();
        s = tmp;
        end_node = s->last_;
    }
    *s = SPath();
    // Rebuild the path.
    *this = SPath(first_);
    assert(last_ == end_node);
    // Update the Node::sp_'s.
    for (Node* n = first_; ; n = n->first_succ()) {
        n->sp_ = this;
        if (n == last_)
            break;
    }
    return true;
}

void Graph::rebuild(const vector<const DNASeq4*>& reads, size_t min_kmer_count) {

    //cerr << "Building graph...\n"; //debug

    clear();

    //
    // Fill the kmer map & count all kmers.
    //
    Kmer km;
    for (const DNASeq4* s : reads) {
        Kmerizer kmers {km_len_, *s};
        Kmer km;
        while ((km = kmers.next()))
            ++map_[km].count;
    }
    //cerr << "Found " << map_.size() << " kmers in " << readset.reads().size() <<" reads.\n"; //debug
    // Remove homopolymers.
    for (Nt2 nt : Nt2::all)
        map_.erase(Kmer::homopolymer(km_len_, nt));

    //
    // Build the standalone nodes.
    //
    for(auto km=map_.begin(); km!=map_.end();) {
        if (km->second.count < min_kmer_count) {
            map_.erase(km++);
        } else {
            nodes_.push_back(Node(NodeData(km->first, km->second.count)));
            // Replace the count of the kmer with the index of the corresponding node.
            km->second.node = nodes_.size() - 1;
            ++km;
        }
    }
    assert(map_.size() == nodes_.size());
    //cerr << "Built " << nodes_.size() << " nodes.\n"; //debug

    //
    // Build the edges.
    //
    for (Node& n : nodes_) {
        assert(!n.empty());
        // Check each possible successor kmer.
        for (size_t nt2=0; nt2<4; ++nt2) {
            auto km = map_.find(n.km().succ(km_len_, nt2));
            if (km != map_.end()) {
                assert(&nodes_[km->second.node] != &n); // Homopolymers were removed.
                assert(km->first.back(km_len_) == Nt2(nt2));
                n.set_succ(nt2, &nodes_[km->second.node]);
            }
        }
    }

    //
    // Build the simple paths.
    //
    for (Node& n : nodes_) {
        assert(!n.empty());
        if (n.n_pred() != 1)
            simple_paths_.push_back(SPath(&n));
        if (n.n_succ() > 1) {
            for (size_t nt2=0; nt2<4; ++nt2) {
                Node* s = n.succ(nt2);
                if (s != NULL && s->n_pred() == 1) {
                    simple_paths_.push_back(SPath(s));
                }
            }
        }
    }
    for (SPath& p : simple_paths_) {
        // A separate loop is required because the vector might have resized.
        assert(!p.empty()); // We've just built the graph; no pieces should have been deleted yet!
        p.update_ptrs();
    }
    //cerr << "Built " << simple_paths_.size() << " simple paths.\n"; //debug
}

bool Graph::remove_cycles() {
    if (!sorted_spaths_.empty())
        // The graph is already acyclic.
        return false;
    return remove_microsat_dimer_cycles();
}

bool Graph::remove_microsat_dimer_cycles() {
    bool edited = false;
    for (SPath& p : simple_paths_) {
        if (p.empty())
            continue;
        if (p.n_nodes() == 1) {
            for (size_t nt2=0; nt2<4; ++nt2) {
                SPath* s = p.succ(nt2);
                if (s == NULL)
                    continue;
                if (s->n_nodes() == 1 && s->succ(size_t(p.last()->km().back(km_len_))) == &p) {
                    if (remove_microsat_dimer_cycle(p, *s))
                        edited = true;
                    break;
                }
            }
        } else if (p.n_nodes() == 2) {
            // This is actually a degenerate case of the above. If no additional paths break the
            // simple path, the two dimer nodes may, depending on the oddness of the kmer
            // size and of the microsatellite tract length, collapse in one single simple
            // path that loops on itself. (e.g. GATATG/k=3: GAT->ATA<->TAT->ATG.)
            for (size_t nt2=0; nt2<4; ++nt2) {
                if (p.succ(nt2) == &p) {
                    p.erase(km_len_);
                    edited = true;
                    break;
                }
            }
        }
    }
    return edited;
}

bool Graph::remove_microsat_dimer_cycle(SPath& p, SPath& q) {
    // Make sure we got this right.
    Nt2 p_front = p.first()->km().front();
    Nt2 q_front = q.first()->km().front();
    Nt2 p_back = p.last()->km().back(km_len_);
    Nt2 q_back = q.last()->km().back(km_len_);
    assert(p.succ(size_t(q_back)) == &q && p.pred(size_t(q_front)) == &q);
    assert(q.succ(size_t(p_back)) == &p && q.pred(size_t(p_front)) == &p);
    // Detect the configuration we're in.
    assert(p.n_succ() >= 1 && p.n_pred() >= 1); // c.f. above.
    assert(q.n_succ() >= 1 && q.n_pred() >= 1);
    size_t p_external = bool(p.n_succ() - 1) + bool(p.n_pred() - 1);
    size_t q_external = bool(q.n_succ() - 1) + bool(q.n_pred() - 1);
    if (p_external == 2 && q_external == 2) {
        // Most general case; we ignore it as there's no clear solution and it's
        // rare anyway.
        #ifdef DEBUG
        cout << "DEBUG: Ignored a microsat dimer that is in the most general configuration.\n";
        #endif
        return false;
    } else if (p_external < 2 && q_external < 2) {
        // Most of these cases correspond to the degenerate case when the two nodes
        // of the cycle are in the same SPath, which is handled as we detect the cycle.
        // DOES_NOT_HAPPEN; -- However, this fails.
        // There's an extra limit case that brings us here -- if the microsat node pair
        // is the extremity of two otherwise disconnected subgraphs going in opposite
        // directions.
        // e.g. AGTGT,CTGT/k=3: AGT->GTG<->TGT<-CTG or TGTGA,TGTC/k=3: TGA<-GTG<->TGT->GTC.
        //
        // In that case, we erase the "first" of the two nodes (whatever that means depends
        // on std::unordered_map); the other one usually will be merged with its in/outcoming
        // simple path).
        assert(( (p.n_succ()-1) == 0 && (q.n_succ()-1) == 0 )
            || ( (p.n_pred()-1) == 0 && (q.n_pred()-1) == 0 ));
        p.erase(km_len_);
        // Note: q.empty() is often true at this point.
    } else if (p_external < 2) {
        assert(q_external == 2);
        p.erase(km_len_);
    } else if (q_external < 2) {
        assert(p_external == 2);
        q.erase(km_len_);
    } else {
        DOES_NOT_HAPPEN;
    }
    return true;
}

bool Graph::topo_sort() const {
    if (!sorted_spaths_.empty())
        // Already sorted.
        return true;
    vector<uchar> visitdata; // 0/1s; whether each spath is a parent in the recursion
    visitdata.reserve(simple_paths_.size());
    for (const SPath& p : simple_paths_)
        if (!p.empty())
            p.visitdata = NULL;

    for (const SPath& p : simple_paths_) {
        if (p.empty())
            continue;
        if(!topo_sort(&p, visitdata)) {
            sorted_spaths_.resize(0);
            return false;
        }
    }
    return true;
}

bool Graph::topo_sort(const SPath* p, vector<uchar>& visitdata) const {
    if (p->visitdata != NULL) {
        if (*(uchar*) p->visitdata)
            // The recursion looped; not a DAG.
            return false;
        else
            // Joining a known path from a different root.
            return true;
    } else {
        visitdata.push_back(true); // n.b. Enough memory was reserved.
        p->visitdata = (void*)&visitdata.back();
        for (size_t nt2=0; nt2<4; ++nt2) {
            const SPath* s = p->succ(nt2);
            if (s != NULL)
                if (!topo_sort(s, visitdata))
                    return false;
        }
        *(uchar*)p->visitdata = false;

        sorted_spaths_.push_back(p);
    }

    return true;
}

bool Graph::find_best_path(vector<const SPath*>& best_path) const {
    best_path.resize(0);

    assert(!empty());
    if(!topo_sort())
        // Not a DAG.
        return false;

    #ifdef DEBUG
    // This is unnecessary as we iterate on a sort.
    for (const SPath& p : simple_paths_)
        if (!p.empty())
            p.visitdata = NULL;
    #endif

    // Compute the best score at each node.
    vector<size_t> scores;
    scores.reserve(sorted_spaths_.size());

    for (const SPath* p : sorted_spaths_) {
        assert(!p->empty());
        //n.b. Terminal nodes were added first.
        size_t succ_scores[4];
        for (size_t nt2=0; nt2<4; ++nt2) {
            const SPath* succ = p->succ(nt2);
            if (succ != NULL)
                //n.b. as the graph is sorted, succ->visitdata has been set.
                succ_scores[nt2] = *(size_t*)succ->visitdata;
            else
                succ_scores[nt2] = 0;
        }
        scores.push_back(*std::max_element(succ_scores, succ_scores+4) + p->km_cumcount());
        p->visitdata = &scores.back();
    }

    // Find the best starting node.
    auto p = sorted_spaths_.rbegin();
    const SPath* best_start = *p;
    size_t best_score = *(size_t*)best_start->visitdata;
    ++p;
    while(p != sorted_spaths_.rend()) {
        if (*(size_t*)(*p)->visitdata > best_score) {
            best_start = *p;
            best_score = *(size_t*)best_start->visitdata;
        }
        ++p;
    }

    // Find the best path.
    const SPath* curr = best_start;
    while(curr != NULL) {
        best_path.push_back(curr);

        size_t succ_scores[4];
        for (size_t nt2=0; nt2<4; ++nt2) {
            const SPath* succ = curr->succ(nt2);
            if (succ != NULL)
                succ_scores[nt2] = *(size_t*)succ->visitdata;
            else
                succ_scores[nt2] = 0;
        }

        size_t* best_succ_score = std::max_element(succ_scores, succ_scores+4);
        curr = curr->succ(best_succ_score - succ_scores);
        // if succ_scores was {0,0,0,0}, curr is null
    }

    return true;
}

void Graph::dump_gfa(const string& path, bool individual_nodes) const {
    ofstream ofs (path);
    if (!ofs) {
        cerr << "Error: Failed to open '" << path << "' for writing.\n";
        throw exception();
    }

    // Write the header.
    ofs << "H\tVN:Z:1.0\n";

    if (!individual_nodes) {
        // Write the vertices.
        for (const SPath& p : simple_paths_)
            if (!p.empty())
                // n.b. In principle the length of the contigs should be (n_nodes+km_len-1).
                // However for visualization purposes we use n_nodes (for now at least).
                ofs << "S\t" << index_of(p.first()) << "\t*\tLN:i:" << p.n_nodes() << "\tKC:i:" << p.km_cumcount() << "\n";

        // Write the edges.
        for (const SPath& p : simple_paths_) {
            if (p.empty())
                continue;
            for (size_t nt2=0; nt2<4; ++nt2) {
                const SPath* succ = p.succ(nt2);
                if (succ != NULL)
                    ofs << "L\t" << index_of(p.first()) << "\t+\t" << index_of(succ->first()) << "\t+\tM" << (km_len_-1) << "\n";
            }
        }
    } else {
        // Write the vertices.
        for (const Node& n : nodes_) {
            if (n.empty())
                continue;
            ofs << "S\t" << index_of(&n) << "\t*\tLN:i:1\tKC:i:" << n.count() << "\tseq:" << n.km().str(km_len_) << "\n";
        }
        // Write the edges.
        for (const Node& n : nodes_) {
            if (n.empty())
                continue;
            for (size_t nt2=0; nt2<4; ++nt2) {
                const Node* succ = n.succ(nt2);
                if (succ != NULL)
                    ofs << "L\t" << index_of(&n) << "\t+\t" << index_of(succ) << "\t+\tM" << (km_len_-1) << "\n";
            }
        }
    }
}

vector<vector<const SPath*>> Graph::components() {
    map<const void*, vector<const SPath*>> components_map;
    for (SPath& p : simple_paths_)
        if (!p.empty())
            p.visitdata = NULL;
    for (SPath& p : simple_paths_) {
        if (p.empty())
            continue;
        propagate_component_id(&p, NULL);
        components_map[p.visitdata].push_back(&p);
    }
    vector<vector<const SPath*>> components;
    for (auto& c : components_map)
        components.push_back(move(c.second));
    return components;
}

void Graph::propagate_component_id(const SPath* p, void* id) {
    if (p->visitdata == NULL) {
        if (id == NULL)
            // (One simple path address serves as an ID for each component.)
            id = const_cast<SPath*>(p);
        p->visitdata = id;
        const SPath* neighbor;
        for (size_t nt2=0; nt2<4; ++nt2) {
            neighbor = p->pred(nt2);
            if (neighbor != NULL)
                propagate_component_id(neighbor, id);
            neighbor = p->succ(nt2);
            if (neighbor != NULL)
                propagate_component_id(neighbor, id);
        }
    }
}
