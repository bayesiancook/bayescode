
#ifndef DIFFSELSPARSEFIT_H
#define DIFFSELSPARSEFIT_H

#include "BidimArray.hpp"

/**
 * \brief An array of site- and condition-specific fitness profiles
 *
 * this array is used in DiffSelSparseModel; It implements the following deterministic relation:
 * - F_0ia = G_0ia;
 * - F_kia = G_0ia^(1-d_kia) * G_kia^(d_kia), for k=1..K-1
 *
 * where
 * - G_kia is the input fitness matrix (for amino-acid a, site i and condition k=0..K-1)
 * - d_kia is an array of toggles (for amino-acid a, site i and condition k=1..K-1)
 * - F_kia is the actual fitness of amino acid a at site i and under condition k (for k=0..K-1)
 *
 * In words, if d_kia == 0, then F_kia, the fitness of amino-acid a at site i under condition k, is just the baseline G_0ia; otherwise, it is a 'new' fitness parameter, such as defined by G_kia.
 * Note that, when Nlevel == 2, the relation between G, d and F unfolds over two levels:
 * - F_0ia = G_0ia;
 * - F_1ia = G_0ia^(1-d_1ia) * G_1ia^(d_1ia);
 * - F_kia = F_1ia^(1-d_kia) * G_kia^(d_kia), for k=2..K-1
 */

class DiffSelSparseFitnessArray : public SimpleBidimArray<vector<double> >    {

    public:

    //! constructor, parameterized by input fitness array, toggle array and Nlevel
    DiffSelSparseFitnessArray(const BidimSelector<vector<double> >& infitness, const BidimSelector<vector<int> >& intoggle, int inNlevel) : 
        SimpleBidimArray<vector<double> >(infitness.GetNrow(),infitness.GetNcol(),vector<double>(infitness.GetVal(0,0).size(),0)),
        fitness(infitness), toggle(intoggle), Nlevel(inNlevel)  {
            Update();
    }

    //! returns dimension of fitness profiles (should normally be 20)
    int GetDim() const {return GetVal(0,0).size();}

    //! full update of the array
    void Update()   {
        for (int i=0; i<GetNrow(); i++)  {
            for (int j=0; j<GetNcol(); j++)   {
                Update(i,j);
            }
        }
    }

    //! update of column j (i.e. site j)
    void UpdateColumn(int j)    {
        for (int i=0; i<GetNrow(); i++) {
            Update(i,j);
        }
    }

    //! update of column j (i.e. site j), but only for those conditions that are flagged
    void UpdateColumn(int j, const vector<int>& flag)  {
        for (int i=0; i<GetNrow(); i++) {
            if (flag[i])    {
                Update(i,j);
            }
        }
    }

    //! update of column j (i.e. site j) and condition i
    void Update(int i, int j)   {
        vector<double>& x = (*this)(i,j);
        double total = 0;
        for (int k=0; k<GetDim(); k++) {
            int l = 0;
            if (i>0)    {
                if (toggle.GetVal(i-1,j)[k])   {
                    l = i;
                }
                else    {
                    if ((Nlevel == 2) && (toggle.GetVal(0,j)[k]))   {
                        l = 1;
                    }
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

