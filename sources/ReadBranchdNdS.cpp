
#include "DistBranchArray.hpp"
#include "dSOmegaPathSuffStat.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char* argv[])    {

    string infile = argv[1];
    int burnin = atoi(argv[2]);
    int until = atoi(argv[3]);
    string treefile = argv[4];
    string outfile = argv[5];

    Tree tree(treefile);
    tree.SetIndices();


    ifstream is(infile.c_str());
    cerr << "burnin\n";
    for (int i=0; i<burnin; i++)    {
        cerr << '.';
        string s;
        getline(is,s);
    }
    cerr << '\n';

    cerr << "reading\n";
    dSOmegaPathSuffStatBranchArray dsomss(tree);
    for (int i=0; i<until-burnin; i++)  {
        cerr << '.';
        dSOmegaPathSuffStatBranchArray tmp(tree);
        is >> tmp;
        dsomss.Add(tmp);
    }
    cerr << '\n';

    DistBranchArray<double> val(tree);
    SimpleBranchArray<double> tmp(tree);
    dsomss.GetdNdS(tmp);
    val.Add(tmp);
    // val.Sort();

    ofstream os((outfile + ".tab").c_str());
    val.TabulateMean(os, 0, 1);

    ofstream tos((outfile + ".tre").c_str());
    val.NodeMeanToStream(tos);

    cerr << "tabulated leaf values in " << outfile << ".tab\n";
    cerr << "newick tree in " << outfile << ".tre\n";
}

