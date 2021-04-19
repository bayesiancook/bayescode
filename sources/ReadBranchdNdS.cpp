
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

void RecursiveBranchToNewick(ostream& os, const Link* from, const BranchSelector<double>& v)    {
    if (from->isLeaf()) {
        os << from->GetNode()->GetName();
    }
    else    {
        os << "(";
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveBranchToNewick(os, link->Out(), v);
            if (link->Next() != from)   {
                os << ",";
            }
        }
        os << ")";
    }
    if (! from->isRoot())	{
        os << ":";
        os << v.GetVal(from->GetBranch()->GetIndex());
    }
}

void BranchToNewick(ostream& os, const BranchSelector<double>& v)   {
    RecursiveBranchToNewick(os, v.GetTree().GetRoot(), v);
    os << ";\n";
}

void RecursiveTabulate(ostream& os, const Link* from, const BranchSelector<double>& v, bool leaf)  {
    if (leaf)   {
        if (from->isLeaf()) {
            os << from->GetNode()->GetName() << '\t' << v.GetVal(from->GetBranch()->GetIndex()) << '\n';
        }
    }
    else    {
        if (! from->isRoot())   {
            os << v.GetTree().GetLeftMost(from) << '\t' << v.GetTree().GetRightMost(from) << '\t' << v.GetVal(from->GetBranch()->GetIndex()) << '\n';
        }
    }
    for (const Link* link=from->Next(); link!=from; link=link->Next())  {
        RecursiveTabulate(os, link->Out(), v, leaf);
    }
}

void Tabulate(ostream& os, const BranchSelector<double>& v, bool leaf) {
    RecursiveTabulate(os, v.GetTree().GetRoot(), v, leaf);
}

int main(int argc, char* argv[])    {

    string suffstatfile = argv[1];
    string treefile = argv[2];
    string outfile = argv[3];

    Tree tree(treefile);
    tree.SetIndices();

    dSOmegaPathSuffStatBranchArray dsomss(tree);
    dSOmegaPathSuffStatBranchArray tmpdsomss(tree);
    ifstream is(suffstatfile.c_str());
    string s;
    int count = 0;
    while (getline(is,s))	{
        istringstream iss(s);
        iss >> tmpdsomss;
        dsomss.Add(tmpdsomss);
        count++;
    }
    cerr << "number of lines read : " << count << '\n';
    SimpleBranchArray<double> dnds(tree);
    dsomss.GetdNdS(dnds);
    SimpleBranchArray<double> ds(tree);
    dsomss.GetdS(ds);

    ofstream tos((outfile + ".dsom.tre").c_str());
    ToNewick(tos, ds, dnds);
    cerr << "newick tree in " << outfile << ".tre\n";

    ofstream tabos((outfile+ ".alldsom.tab").c_str());
    Tabulate(tabos, dnds, false);
    ofstream ltabos((outfile+ ".leafdsom.tab").c_str());
    Tabulate(ltabos, dnds, true);
    cerr << "tabulated branch dN/dS values in : " << outfile << ".dsom.tab\n";
    cerr << "for terminal branches only       : " << outfile << ".leafdsom.tab\n";

    SimpleBranchArray<double> dscount(tree);
    dsomss.GetSynCount(dscount);
    SimpleBranchArray<double> dsbeta(tree);
    dsomss.GetSynBeta(dsbeta);
    SimpleBranchArray<double> dncount(tree);
    dsomss.GetNonSynCount(dncount);
    SimpleBranchArray<double> dnbeta(tree);
    dsomss.GetNonSynBeta(dnbeta);

    ofstream dScos((outfile + ".counts_dS.dnd").c_str());
    BranchToNewick(dScos, dscount);
    ofstream dSbos((outfile + ".counts_dS_norm.dnd").c_str());
    BranchToNewick(dSbos, dsbeta);
    ofstream dNcos((outfile + ".counts_dN.dnd").c_str());
    BranchToNewick(dNcos, dncount);
    ofstream dNbos((outfile + ".counts_dN_norm.dnd").c_str());
    BranchToNewick(dNbos, dnbeta);
    cerr << "suffstats in " << outfile << ".counts_dX[_norm].dnd\n";
}

