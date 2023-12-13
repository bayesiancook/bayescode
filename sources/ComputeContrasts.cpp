#include "dSOmegaPathSuffStat.hpp"
#include "RecursiveNewick.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char* argv[])    {

    // continuous data
    ifstream is(argv[1]);
    Tree tree(argv[2]);
    string basename = argv[3];

    tree.SetIndices();

    dSOmegaPathSuffStatBranchArray dsomss(tree);

    int ngene;
    is >> ngene;
    for (int gene=0; gene<ngene; gene++)    {
        string genename;
        is >> genename;

        string tmp;
        is >> tmp;
        if (tmp != "counts_dS") {
            cerr << "error when reading suffstat file\n";
            exit(1);
        }
        Tree treedscount(is);
        treedscount.SetIndices();

        is >> tmp;
        if (tmp != "counts_dS_norm") {
            cerr << "error when reading suffstat file\n";
            exit(1);
        }
        Tree treedsbeta(is);
        treedsbeta.SetIndices();

        is >> tmp;
        if (tmp != "counts_dN") {
            cerr << "error when reading suffstat file\n";
            exit(1);
        }
        Tree treedncount(is);
        treedncount.SetIndices();

        is >> tmp;
        if (tmp != "counts_dN_norm") {
            cerr << "error when reading suffstat file\n";
            exit(1);
        }
        Tree treednbeta(is);
        treednbeta.SetIndices();

        dsomss.Add(treedscount, treedsbeta, treedncount, treednbeta);
    }

    SimpleBranchArray<double> dnds(tree);
    dsomss.GetdNdS(dnds);
    SimpleBranchArray<double> ds(tree);
    dsomss.GetdS(ds);

    ofstream dsom_os((basename + ".dsompathsuffstat").c_str());
    dsom_os << "counts_dS\t";
    dsomss.BranchToNewickSynCount(dsom_os);
    dsom_os << "counts_dS_norm\t";
    dsomss.BranchToNewickSynBeta(dsom_os);
    dsom_os << "counts_dN\t";
    dsomss.BranchToNewickNonSynCount(dsom_os);
    dsom_os << "counts_dN_norm\t";
    dsomss.BranchToNewickNonSynBeta(dsom_os);

    ofstream tos((basename + ".dsom.tre").c_str());
    ToNewick(tos, ds, dnds);
    cerr << "newick tree in " << basename << ".dsom.tre\n";

    ofstream tabos((basename + ".alldsom.tab").c_str());
    Tabulate(tabos, ds, dnds, false);
    ofstream ltabos((basename + ".leafdsom.tab").c_str());
    Tabulate(ltabos, ds, dnds, true);
    cerr << "tabulated branch dS and dN/dS in : " << basename << ".dsom.tab\n";
    cerr << "for terminal branches only       : " << basename << ".leafdsom.tab\n";

}
