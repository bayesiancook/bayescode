#ifndef TAXONSET_H
#define TAXONSET_H

#include <iostream>
#include <map>
#include <string>
#include <vector>
#define nullptr 0 

class Tree;
class Link;

class TaxonSet {
  public:
    TaxonSet(const std::vector<std::string>& names);
    TaxonSet(const Tree* tree, const Link *subgroup = nullptr);
    ~TaxonSet() = default;

    int GetNtaxa() const {return Ntaxa;}
    std::string GetTaxon(int index) const {return taxlist[index];}
    int GetTaxonIndex(std::string intaxon) const;
    int GetTaxonIndexWithIncompleteName(std::string taxname) const;

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
