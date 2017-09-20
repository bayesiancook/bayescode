
#ifndef DIFFSELFIT_H
#define DIFFSELFIT_H

#include "BidimArray.hpp"

class DiffSelFitnessArray : public SimpleBidimArray<vector<double> >    {

    public:

    DiffSelFitnessArray(const ConstArray<vector<double> >& inbaseline, const ConstBidimArray<vector<double> >& indelta, int inNlevel) : 
    // DiffSelFitnessArray(const ConstArray<vector<double> >& inbaseline, const ConstBidimArray<vector<double> >& indelta, const vector<vector<int> >& inpattern) :
        SimpleBidimArray<vector<double> >(indelta.GetNrow()+1,indelta.GetNcol(),vector<double>(indelta.GetVal(0,0).size(),0)),
        baseline(inbaseline), delta(indelta), Nlevel(inNlevel)  {
        // baseline(inbaseline), delta(indelta), pattern(inpattern) {
            Update();
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
            double d = 0;
            if ((i==1) || (Nlevel == 2))   {
                d += delta.GetVal(0,j)[k];
            }
            if (i > 1)  {
                d += delta.GetVal(i-1,j)[k];
            }
            /*
            for (unsigned int l=1; l<pattern[i].size(); l++)    {
                if (pattern[l][i])  {
                    d += delta.GetVal(l-1,j)[k];
                }
            }
            */
            x[k] = baseline.GetVal(j)[k] * exp(d);
            total += x[k];
        }
        for (int k=0; k<GetDim(); k++) {
            x[k] /= total;
        }
    }

    protected:

    const ConstArray<vector<double> >& baseline;
    const ConstBidimArray<vector<double> >& delta;
    // const vector<vector<int> >& pattern;
    int Nlevel;
};

#endif

