#ifndef TAXONSET_H
#define TAXONSET_H

#include <iostream>
#include <map>
#include <string>
#include <vector>
#define nullptr 0 

class Tree;
class Link;

/**
 * \brief A set of taxa (simply, a bi-directional correspondance between taxon names and integers between 0 and Ntaxa-1)
 *
 * TaxonSet implements a globally consistent correspondance of taxon names and indices between trees and sequence alignments.
 * Each Tree and SequenceAlignment object has a pointer to a TaxonSet.
 * Typically, a SequenceAlignment is created by reading from a file;
 * the TaxonSet is created at that step, by the SequenceAlignment constructor itself.
 * Then, a Tree is read from file and is registered with the taxon set given by the sequence alignment.
 * Once this is done, both objects have mutually consistent indexing systems and can work together (e.g. to calculate a likelihood).
 *
 * Of note: some work to do here to ensure proper destruction of TaxonSet (smart pointers?)
 */

class TaxonSet {
  public:

    //! constructor, based on a vector of taxon names
    TaxonSet(const std::vector<std::string>& names);
    //! constructor, based on a vector of taxon names
    TaxonSet(const TaxonSet& from);
    //! default constructor
    ~TaxonSet() = default;

    //! return number of taxa
    int GetNtaxa() const {return Ntaxa;}
    //! return taxon name, given the index
    std::string GetTaxon(int index) const {return taxlist[index];}
    //! return taxon index, given the name
    int GetTaxonIndex(std::string intaxon) const;
    //! return taxon index, given incomplete name (first part)
    int GetTaxonIndexWithIncompleteName(std::string taxname) const;
    //! formatted output to stream
    void ToStream(std::ostream &os) const;

  private:
    void RecursiveGetSubSet(const Link *from, int &i);

    int Ntaxa;
    std::map<std::string, int> taxmap;
    std::vector<std::string> taxlist;
};

inline int TaxonSet::GetTaxonIndex(std::string intaxon) const { 
    std::map<std::string,int>::const_iterator i = taxmap.find(intaxon);
    if (i == taxmap.end())  {
        std::cerr << "error in TaxonSet: taxon not found\n";
        exit(1);
    }
    return i->second - 1; 
}

#endif  // TAXONSET_H
