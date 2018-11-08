#pragma once

#include <set>
#include "CodonSequenceAlignment.hpp"

/**
 * \brief A data wrapper around the polymorphism data, created from a codon sequence
 * alignment and an absolute path to a alignment.
 * PolyData will search automatically for polymorphism data (.vcf files)
 * that are in the same directory (folder) as the alignment.
 * The data will be match to the taxon (leaf of the tree) based on file names,
 * meaning the file (.vcf) must include the taxon name.
 */

class PolyData {
  public:
    //! \brief Constructor: takes a (const pointer to a) codon sequence
    //! alignment and an absolute path to a file (string).
    //! PolyData will search for every .vcf files that are in the same
    //! directory (folder) as ali_path. The .vcf files  will be matched
    //! to taxa (leaves of the tree) if the name of the taxon is included
    //! in the file name (e.g. *HomoSapiens*.vcf)
    PolyData(CodonSequenceAlignment *from, std::string const &ali_path);

    ~PolyData() /*override*/ = default;

    //! return size of state space
    int GetNstate() const;

    //! return the number of taxa (number of aligned sequences)
    int GetNtaxa() const;

    //! return the number of aligned positions
    int GetNsite() const;

    //! return the number of files (.vcf) found
    int GetNvcf() const;

    //! return the number of occurrence of allele (state) in leaf (taxon)
    //! at position (site) of the alignment, return 0 if missing entry
    //! return 1 is the position is not polymorphic
    unsigned GetCount(int taxon, int site, int state) const;

    //! return the number of sampled allele in leaf (taxon)
    //! at position (site) of the alignment, return 0 if missing entry,
    //! return 1 is the position is not polymorphic
    unsigned GetSampleSize(int taxon, int site) const;

    std::set<unsigned> GetSampleSizeSet() const {
        return std::set<unsigned>(SampleSize.begin(), SampleSize.end());
    };

  private:
    //! a sparse data collection for the count (number of copies) of alleles,
    //! for each taxa (1st level), for each site of the alignment (2nd level),
    //! and for each state of the codons (3rd level).
    std::map<int, std::map<int, std::map<int, unsigned>>> Data;

    //! a sparse data collection for the total count of sampled allele,
    //! for each taxa (1st level), and for each site of the alignment (2nd level).
    std::map<int, std::map<int, unsigned>> SiteSampleSize;

    //! a map for the total count of sampled allele for each taxa.
    std::vector<unsigned> SampleSize;

    //! a codon sequence alignment
    CodonSequenceAlignment *Alignment;

    //! The number of files (.vcf) found
    int Nvcf{0};
};
