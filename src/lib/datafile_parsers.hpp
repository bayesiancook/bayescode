#pragma once

#include <fstream>
#include <map>
#include <set>
#include <string>
#include "SequenceAlignment.hpp"

using GeneSet = std::vector<std::string>;
using GeneLengths = std::map<std::string, int>;

GeneSet parse_datafile(std::string filename) {
    std::ifstream is(filename);
    int nb_genes;
    is >> nb_genes;
    GeneSet result;
    for (int i = 0; i < nb_genes; i++) {
        std::string name;
        is >> name;
        result.push_back(name);
    }
    return result;
}

GeneLengths parse_geneset_alignments(const GeneSet& gene_set) {
    GeneLengths result;
    for (auto gene : gene_set) {
        FileSequenceAlignment tmp(gene);
        result[gene] = tmp.GetNsite() / 3;
    }
    return result;
}