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

