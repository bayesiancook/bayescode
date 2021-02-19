#pragma once
#include "BranchArray.hpp"
#include <list>

template<class T> class DistBranchArray   {

    public:

    DistBranchArray(const Tree& intree) : tree(intree), branchdist(intree.GetNbranch(), list<T>()), counter(0)  {}

    virtual const Tree &GetTree() const { return tree; }
    const Link* GetRoot() const { return tree.GetRoot(); }
    int GetNbranch() const { return tree.GetNbranch(); }

    void Add(const BranchSelector<T>& datapoint)    {
        for (int i=0; i<GetNbranch(); i++) {
            branchdist[i].push_front(datapoint.GetVal(i));
        }
        counter++;
    }

    void Sort() {
        for (int i=0; i<GetNbranch(); i++) {
            branchdist[i].sort();
        }
    }

    T GetMean(int index) const  {
		const list<T>& l = branchdist[index];
        T mean = 0;
		for (auto val : l) {
            mean += val; 
		}
        mean /= l.size();
        return mean;
    }

    const T& GetQuantile(int index, double q) const {
		const list<T>& l = branchdist[index];
		auto i = l.begin();
		int n = ((int) (((double) l.size())*q));
		for (int j=0; j<n; j++)	{
			i++;
		}
        return *i;
    }

    void TabulateMean(ostream& os) const	{
        os << "#NodeName\tmean\tmedian\tmin95\tmax95\n";
        RecursiveTabulateMean(os, GetRoot());
    }

    void RecursiveTabulateMean(ostream& os, const Link* from) const {
        if (! from->isRoot())   {
            int index = from->GetBranch()->GetIndex();
            os << from->GetNode()->GetName() << '\t' << GetMean(index) << '\t' << GetQuantile(index, 0.5) << '\t' << GetQuantile(index, 0.025) << '\t' << GetQuantile(index, 0.975) << '\n';
        }
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveTabulateMean(os, link->Out());
        }
    }

    void MeanToStream(ostream& os) const {
        RecursiveMeanToStream(os, GetRoot());
        os << ";\n";
    }

    void RecursiveMeanToStream(ostream& os, const Link* from) const {
        if (from->isLeaf()) {
            os << from->GetNode()->GetName();
        }
        else    {
            os << "(";
            for (const Link* link=from->Next(); link!=from; link=link->Next())  {
                RecursiveMeanToStream(os,link->Out());
                if (link->Next() != from)   {
                    os << ",";
                }
            }
            os << ")";
	    os << from->GetNode()->GetName();
        }
        if (! from->isRoot())    {
            os << ":";
            os << GetMean(from->GetBranch()->GetIndex());
        }
    }

    private:

    const Tree& tree;
    vector<list<T>> branchdist;
    int counter;
};

