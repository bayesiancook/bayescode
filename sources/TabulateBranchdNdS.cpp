#include "dSOmegaPathSuffStat.hpp"
#include "RecursiveNewick.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char* argv[])    {

    ifstream is(argv[1]);
    Tree tree(argv[2]);
    tree.SetIndices();

    string prefix = argv[3];
    string basename = argv[4];

    dSOmegaPathSuffStatBranchArray dsomss(tree);

    int ngene;
    is >> ngene;
    for (int gene=0; gene<ngene; gene++)    {
        string genename;
        is >> genename;

        Tree treedscount(prefix + genename + ".counts_dS.dnd");
        treedscount.SetIndices();

        Tree treedsbeta(prefix + genename + ".counts_dS_norm.dnd");
        treedsbeta.SetIndices();

        Tree treedncount(prefix + genename + ".counts_dN.dnd");
        treedncount.SetIndices();

        Tree treednbeta(prefix + genename + ".counts_dN_norm.dnd");
        treednbeta.SetIndices();

        dsomss.Add(treedscount, treedsbeta, treedncount, treednbeta);
    }

    ofstream dos((basename + ".branchdsomsuffstat").c_str());
    dos << dsomss;
    cerr << "tabulated suff stats in " << basename << ".branchdsomsuffstat\n";
    
    SimpleBranchArray<double> dnds(tree);
    dsomss.GetdNdS(dnds);
    SimpleBranchArray<double> ds(tree);
    dsomss.GetdS(ds);

    ofstream tos((basename + ".dsom.tre").c_str());
    ToNewick(tos, ds, dnds);
    cerr << "newick tree in " << basename << ".dsom.tre\n";

    ofstream tabos((basename + ".alldsom.tab").c_str());
    Tabulate(tabos, dnds, false);
    ofstream ltabos((basename + ".leafdsom.tab").c_str());
    Tabulate(ltabos, dnds, true);
    cerr << "tabulated branch dN/dS values in : " << basename << ".dsom.tab\n";
    cerr << "for terminal branches only       : " << basename << ".leafdsom.tab\n";

}

