#pragma once
#include "Chronogram.hpp"

template<class T> class DistChronoBranchArray   {

    public:

    DistChronoBranchArray(const Tree& intree) : tree(intree), 
        branchmean(intree.GetNbranch(), 0), 
        branchtime(intree.GetNbranch(), 0), counter(0)  {
    }

    virtual const Tree &GetTree() const { return tree; }
    const Link* GetRoot() const { return tree.GetRoot(); }
    int GetNbranch() const { return tree.GetNbranch(); }

    void Add(const Chronogram& chrono, const BranchSelector<T>& datapoint)    {
        RecursiveAdd(GetRoot(), chrono, datapoint);
        counter++;
    }

    void RecursiveAdd(const Link* from, const Chronogram& chrono, const BranchSelector<T>& datapoint)    {
        if (! from->isRoot())   {
            branchmean[from->GetBranch()->GetIndex()] += datapoint.GetVal(from->GetBranch()->GetIndex());
            double dt = chrono.GetDeltaTime(from);
            branchtime[from->GetBranch()->GetIndex()] += dt;
        }
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveAdd(link->Out(), chrono, datapoint);
        }
    }

    void Normalize()    {
        for (int i=0; i<GetNbranch(); i++) {
            branchmean[i] /= counter;
            branchtime[i] /= counter;
        }
    }

    void Tabulate(ostream& os) const    {
        for (int i=0; i<GetNbranch(); i++) {
            os << branchtime[i] << '\t' << branchmean[i] << '\n';
        }
    }

    private:

    const Tree& tree;
    vector<T> branchmean;
    vector<double> branchtime;
    int counter;
};

