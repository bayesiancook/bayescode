#include <cstdlib>
#include <iostream>
using namespace std;

#include "BiologicalSequences.hpp"
#include "TaxonSet.hpp"
#include "tree/implem.hpp"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//     TaxonSet
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

TaxonSet::TaxonSet(const std::vector<string> &names) : Ntaxa(names.size()), taxlist(names) {
    for (int i = 0; i < Ntaxa; i++) {
        if (taxmap[names[i]] != 0) {
            cerr << "found several taxa with same name : " << names[i] << '\n';
            exit(1);
        }
        taxmap[names[i]] = i + 1;
    }
}

std::vector<int> TaxonSet::get_reverse_index_table(const Tree *tree) const {
    std::vector<int> ret(0, -1);
    for (size_t node = 0; node < tree->nb_nodes(); node++) {
        if (tree->is_leaf(node)) {
            ret.emplace(ret.begin() + GetTaxonIndex(tree->node_name(node)), node);
        }
    }
    return ret;
}

std::vector<int> TaxonSet::get_index_table(const Tree *tree) const {
    std::vector<int> ret(tree->nb_nodes(), -1);
    for (size_t node = 0; node < tree->nb_nodes(); node++) {
        if (tree->is_leaf(node)) {
            if (tree->node_name(node) == "") {
                cerr << "error: leaf has no name\n";
                exit(1);
            }
            ret[node] = GetTaxonIndex(tree->node_name(node));
        }
    }
    cerr << "get index table ok\n";
    return ret;
}

TaxonSet::TaxonSet(const TaxonSet &from)
    : Ntaxa(from.GetNtaxa()), taxmap(from.taxmap), taxlist(from.taxlist) {}

void TaxonSet::ToStream(ostream &os) const {
    os << Ntaxa << '\n';
    for (int i = 0; i < Ntaxa; i++) { os << taxlist[i] << '\n'; }
}

int TaxonSet::GetTaxonIndexWithIncompleteName(string taxname) const {
    int found = -1;
    for (int i = 0; i < Ntaxa; i++) {
        if (taxlist[i].substr(0, taxname.length()) == taxname) {
            if (found != -1) {
                cerr << "error : taxon found twice : " << taxname << '\n';
                exit(1);
            }
            found = i;
        }
    }
    if (found == -1) {
        for (int i = 0; i < Ntaxa; i++) {
            if (taxname.substr(0, taxlist[i].length()) == taxlist[i]) {
                if (found != -1) {
                    cerr << "error : taxon found twice : " << taxname << '\n';
                    exit(1);
                }
                found = i;
            }
        }
    }
    return found;
}
