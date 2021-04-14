
#include "DistBranchArray.hpp"
#include "dSOmegaPathSuffStat.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

void RecursiveToNewick(ostream& os, const Link* from, const BranchSelector<double>& ds, const BranchSelector<double>& dnds) {
    if (from->isLeaf()) {
        os << from->GetNode()->GetName();
        os << "_";
    }
    else    {
        os << "(";
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveToNewick(os, link->Out(), ds, dnds);
            if (link->Next() != from)   {
                os << ",";
            }
        }
        os << ")";
    }
    if (! from->isRoot())	{
        os << dnds.GetVal(from->GetBranch()->GetIndex());
        os << ":";
        os << ds.GetVal(from->GetBranch()->GetIndex());
    }
    else    {
        os << dnds.GetVal(from->Next()->GetBranch()->GetIndex());
    }
}

void ToNewick(ostream& os, const BranchSelector<double>& ds, const BranchSelector<double>& dnds)    {
    RecursiveToNewick(os, ds.GetTree().GetRoot(), ds, dnds);
    os << ";\n";
}

int main(int argc, char* argv[])    {

    string suffstatfile = argv[1];
    string treefile = argv[2];
    string outfile = argv[3];

    Tree tree(treefile);
    tree.SetIndices();

    dSOmegaPathSuffStatBranchArray dsomss(tree);
    ifstream is(suffstatfile.c_str());
    is >> dsomss;
    SimpleBranchArray<double> dnds(tree);
    dsomss.GetdNdS(dnds);
    SimpleBranchArray<double> ds(tree);
    dsomss.GetdS(ds);

    ofstream tos((outfile + ".tre").c_str());
    ToNewick(tos, ds, dnds);
    cerr << "newick tree in " << outfile << ".tre\n";

    /*
    ofstream os((outfile + ".tab").c_str());
    val.TabulateMean(os, 0, 1);
    cerr << "tabulated leaf values in " << outfile << ".tab\n";
    */
}

