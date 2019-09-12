#pragma once
#include <fstream>
#include "global/logging.hpp"
#include "lib/CodonSequenceAlignment.hpp"
#include "tree/implem.hpp"

struct PreparedData {
    NHXParser parser;
    std::unique_ptr<const Tree> tree;
    FileSequenceAlignment nuc_align;
    CodonSequenceAlignment alignment;
    TaxonSet taxon_set;

    PreparedData(std::string alignfile, std::ifstream& treefile)
        : parser(treefile),
          tree(make_from_parser(parser)),
          nuc_align(alignfile),
          alignment(&nuc_align),
          taxon_set(*alignment.GetTaxonSet()) {
        // various checks and debug
        assert(tree->nb_nodes() > 0);
        DEBUG("Parsed tree with {} nodes.", tree->nb_nodes());

        assert(nuc_align.GetNtaxa() > 0 && nuc_align.GetNsite() > 0);
        DEBUG("Parsed alignment with {} sequences of length {}. Example taxon name: {}.",
            nuc_align.GetNtaxa(), nuc_align.GetNsite(), nuc_align.GetTaxonSet()->GetTaxon(0));

        assert(alignment.GetNtaxa() > 0 && alignment.GetNsite() > 0);
        DEBUG("Converted alignment to codons (new length: {}).", alignment.GetNsite());

        DEBUG("Got a taxon set of length {}. Example taxon name: {}.", taxon_set.GetNtaxa(),
            taxon_set.GetTaxon(0));
    }
};

PreparedData prepare_data(std::string alignfile, std::string treefile) {
    // parsing tree
    std::ifstream tree_stream{treefile};
    return {alignfile, tree_stream};
}