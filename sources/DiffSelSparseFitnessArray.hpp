
#ifndef DIFFSELSPARSEFIT_H
#define DIFFSELSPARSEFIT_H

#include "BidimArray.hpp"

class DiffSelSparseFitnessArray : public SimpleBidimArray<vector<double> >    {

    public:

    DiffSelSparseFitnessArray(const BidimSelector<vector<double> >& infitness, const BidimSelector<vector<int> >& intoggle, int inNlevel) : 
        SimpleBidimArray<vector<double> >(infitness.GetNrow(),infitness.GetNcol(),vector<double>(infitness.GetVal(0,0).size(),0)),
        fitness(infitness), toggle(intoggle), Nlevel(inNlevel)  {
            cerr << "dim of fitness array: " << GetNrow() << '\t' << GetNcol() << '\t' << GetDim() << '\n';
            Update();
            cerr << "update ok\n";
    }

    int GetDim() const {return GetVal(0,0).size();}

    void Update()   {
        for (int i=0; i<GetNrow(); i++)  {
            for (int j=0; j<GetNcol(); j++)   {
                Update(i,j);
            }
        }
    }

    void UpdateColumn(int j)    {
        for (int i=0; i<GetNrow(); i++) {
            Update(i,j);
        }
    }

    void UpdateColumn(int j, const vector<int>& flag)  {
        for (int i=0; i<GetNrow(); i++) {
            if (flag[i])    {
                Update(i,j);
            }
        }
    }

    void Update(int i, int j)   {
        vector<double>& x = (*this)(i,j);
        double total = 0;
        for (int k=0; k<GetDim(); k++) {
            int l = 0;
            if (i>0)    {
                if ((Nlevel == 2) && (toggle.GetVal(0,j)[k]))   {
                    l = 1;
                }
                if ((i>1) && (toggle.GetVal(i-1,j)[k]))   {
                    l = i;
                }
            }
            x[k] = fitness.GetVal(l,j)[k];
            total += x[k];
        }
        for (int k=0; k<GetDim(); k++) {
            x[k] /= total;
        }
    }

    protected:

    const BidimSelector<vector<double> >& fitness;
    const BidimSelector<vector<int> >& toggle;
    int Nlevel;
};

#endif

