#pragma once

#include "BranchArray.hpp"
#include "NodeArray.hpp"
#include <list>

class DistBranchNodeArray   {

    public:

    DistBranchNodeArray(const Tree& intree) : tree(intree), branchdist(intree.GetNbranch(), list<double>()), nodedist(intree.GetNnode(), list<double>()), counter(0)  {}

    virtual const Tree &GetTree() const { return tree; }
    const Link* GetRoot() const { return tree.GetRoot(); }
    int GetNbranch() const { return tree.GetNbranch(); }
    int GetNnode() const { return tree.GetNnode(); }

    void Add(const BranchSelector<double>& branchval, const NodeSelector<vector<double>>& nodeval, int index)    {
        for (int i=0; i<GetNbranch(); i++) {
            branchdist[i].push_front(branchval.GetVal(i));
        }
        for (int i=0; i<GetNnode(); i++)    {
            nodedist[i].push_front(nodeval.GetVal(i)[index]);
        }
        counter++;
    }

    void AddFromChrono(const NodeSelector<double>& chrono, const NodeSelector<vector<double>>& nodeval, int index)    {
        RecursiveAddFromChrono(GetRoot(), chrono, nodeval, index);
        counter++;
    }

    void RecursiveAddFromChrono(const Link* from, const NodeSelector<double>& chrono, const NodeSelector<vector<double>>& nodeval, int index)    {
        if (! from->isRoot())   {
            double bl = chrono.GetVal(from->Out()->GetNode()->GetIndex()) - chrono.GetVal(from->GetNode()->GetIndex());
            branchdist[from->GetBranch()->GetIndex()].push_front(bl);
        }
        nodedist[from->GetNode()->GetIndex()].push_front(nodeval.GetVal(from->GetNode()->GetIndex())[index]);

        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            RecursiveAddFromChrono(link->Out(), chrono, nodeval, index);
        }
    }

    void Sort() {
        for (int i=0; i<GetNbranch(); i++) {
            branchdist[i].sort();
        }
        for (int i=0; i<GetNnode(); i++)    {
            nodedist[i].sort();
        }
    }

    double GetBranchMean(int index) const  {
		const list<double>& l = branchdist[index];
        double mean = 0;
		for (auto val : l) {
            mean += val; 
		}
        mean /= l.size();
        return mean;
    }

    double GetBranchQuantile(int index, double q) const {
		const list<double>& l = branchdist[index];
		auto i = l.begin();
		int n = ((int) (((double) l.size())*q));
		for (int j=0; j<n; j++)	{
			i++;
		}
        return *i;
    }

    double GetNodeMean(int index) const  {
		const list<double>& l = nodedist[index];
        double mean = 0;
		for (auto val : l) {
            mean += val; 
		}
        mean /= l.size();
        return mean;
    }

    double GetNodeQuantile(int index, double q) const {
		const list<double>& l = nodedist[index];
		auto i = l.begin();
		int n = ((int) (((double) l.size())*q));
		for (int j=0; j<n; j++)	{
			i++;
		}
        return *i;
    }

    void MedianToStream(ostream& os) const {
        RecursiveMedianToStream(os, GetRoot());
        os << ";\n";
    }

    void RecursiveMedianToStream(ostream& os, const Link* from) const {
        if (from->isLeaf()) {
            os << from->GetNode()->GetName();
            os << "_";
        }
        else    {
            os << "(";
            for (const Link* link=from->Next(); link!=from; link=link->Next())  {
                RecursiveMedianToStream(os,link->Out());
                if (link->Next() != from)   {
                    os << ",";
                }
            }
            os << ")";
        }
        os << exp(GetNodeQuantile(from->GetNode()->GetIndex(), 0.5));
        if (! from->isRoot())    {
            os << ":";
            os << GetBranchMean(from->GetBranch()->GetIndex());
        }
    }

    private:

    const Tree& tree;
    vector<list<double>> branchdist;
    vector<list<double>> nodedist;
    int counter;
};

