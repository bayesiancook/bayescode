
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
    string genefile = argv[3];
    double cuts = atof(argv[4]);
    double cutn = atof(argv[5]);
    string outfile = argv[6];

    Tree tree(treefile);
    tree.SetIndices();

    ifstream gis(genefile.c_str());
    int Ngene;
    gis >> Ngene;
    vector<string> genename(Ngene);
    for (int i=0; i<Ngene; i++) {
        int tmp;
        gis >> genename[i] >> tmp;
    }

    vector<dSOmegaPathSuffStatBranchArray> genedsomss(Ngene, dSOmegaPathSuffStatBranchArray(tree));
    ifstream is(suffstatfile.c_str());
    for (int i=0; i<Ngene; i++) {
        dSOmegaPathSuffStatBranchArray tmp(tree);
        is >> tmp;
        genedsomss[i].Add(tmp);
    }

    // get gene omega's
    vector<double> omega(Ngene, 0);
    for (int i=0; i<Ngene; i++) {
        omega[i] = genedsomss[i].GetMeandNdS();
    }

    // get branch modulators
    dSOmegaPathSuffStatBranchArray weighteddsomss(tree);
    for (int i=0; i<Ngene; i++) {
        weighteddsomss.Add(genedsomss[i], omega[i]);
    }

    SimpleBranchArray<double> beta(tree);
    weighteddsomss.GetdNdS(beta);

    SimpleBranchArray<double> ds(tree);
    weighteddsomss.GetdS(ds);

    ofstream zos((outfile + ".z").c_str());
    dSOmegaPathSuffStatBranchArray totdsomss(tree);
    int nincluded = 0;
    for (int i=0; i<Ngene; i++) {
        // SimpleBranchArray<double> zs(tree);
        double maxzs = genedsomss[i].GetMaxZS(ds);
        double maxzn = genedsomss[i].GetMaxZN(ds, beta, omega[i]);
        zos << i << '\t' << genename[i] << '\t' << omega[i] << '\t' << maxzs << '\t' << maxzn << '\n';
        if ((maxzs < cuts) && (maxzn < cutn))   {
            nincluded ++;
            totdsomss.Add(genedsomss[i]);
        }
    }

    cerr << nincluded << " genes selected out of " << Ngene << '\n';

    SimpleBranchArray<double> sel_dnds(tree);
    totdsomss.GetdNdS(sel_dnds);

    SimpleBranchArray<double> sel_ds(tree);
    totdsomss.GetdS(sel_ds);

    ofstream tos((outfile + ".dsom.tre").c_str());
    ToNewick(tos, sel_ds, sel_dnds);
    cerr << "newick tree in " << outfile << ".dsom.tre\n";

    ofstream tabos((outfile+ ".alldsom.tab").c_str());
    Tabulate(tabos, sel_dnds, false);
    ofstream ltabos((outfile+ ".leafdsom.tab").c_str());
    Tabulate(ltabos, sel_dnds, true);
    cerr << "tabulated branch dN/dS values in : " << outfile << ".alldsom.tab\n";
    cerr << "for terminal branches only       : " << outfile << ".leafdsom.tab\n";

    ofstream ssos((outfile + ".branchdsomsuffstat").c_str());
    ssos << totdsomss << '\n';

    /*
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
    */
}

