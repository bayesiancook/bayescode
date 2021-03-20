
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

    tree.ToStream(cerr);
    cerr << '\n';

    cerr << "create\n";
    DistBranchArray<double> val(tree);

    ifstream is(infile.c_str());
    cerr << "burnin\n";
    for (int i=0; i<burnin; i++)    {
        cerr << '.';
        string s;
        getline(is,s);
    }
    cerr << '\n';

    cerr << "reading\n";
    for (int i=0; i<until-burnin; i++)  {
        cerr << '.';
        dSOmegaPathSuffStatBranchArray dsomss(tree);
        is >> dsomss;
        SimpleBranchArray<double> tmp(tree);
        dsomss.GetdNdS(tmp);
        val.Add(tmp);
    }
    cerr << '\n';

    cerr << "sort\n";
    val.Sort();

    ofstream os(outfile.c_str());
    cerr << "tabulate\n";
    val.TabulateMean(os, 0, 1);
    cerr << "ok\n";
}

