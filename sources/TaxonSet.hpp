#ifndef TAXONSET_H
#define TAXONSET_H

#include <iostream>
#include <map>
#include <string>

class Tree;
class Link;

class TaxonSet {
  public:
    TaxonSet(const std::string *names, int ntaxa);
    TaxonSet(const Tree *tree, const Link *subgroup = nullptr);
    ~TaxonSet();

    int GetNtaxa() const;
    std::string GetTaxon(int index) const;
    int GetTaxonIndex(std::string intaxon) const;
    int GetTaxonIndexWithIncompleteName(std::string taxname) const;

    void ToStream(std::ostream &os) const;

  private:
    void RecursiveGetSubSet(const Link *from, int &i);

    int Ntaxa;
    std::map<std::string, int> taxmap;
    std::string *taxlist;
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

inline int TaxonSet::GetNtaxa() const { return Ntaxa; }
inline std::string TaxonSet::GetTaxon(int index) const { return taxlist[index]; }
inline int TaxonSet::GetTaxonIndex(std::string intaxon) const { 
    std::map<std::string,int>::const_iterator i = taxmap.find(intaxon);
    if (i == taxmap.end())  {
        std::cerr << "error in TaxonSet: taxon not found\n";
        exit(1);
    }
    return i->second - 1; 
}

#endif  // TAXONSET_H
