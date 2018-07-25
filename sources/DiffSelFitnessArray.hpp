
#ifndef DIFFSELFIT_H
#define DIFFSELFIT_H

#include "BidimArray.hpp"

/**
 * \brief An array of site- and condition-specific fitness profiles
 *
 * this array is used in DiffSelModel; It implements the following deterministic
 * relation:
 * - F_0ia = G_ia;
 * - F_kia = G_ia * exp(D_kia), for k=1..K-1
 *
 * where G_ia is the baseline fitness of amino-acid a at site i, D_kia the
 * differential effect for amino-acid a at site i under condition k; the result,
 * F_kia is the fitness of amino-acid a at site i under condition k. Note that,
 * when Nlevel == 2, the relation between G, D and F unfolds over two levels:
 * - F_0ia = G_ia;
 * - F_1ia = G_ia * exp(D_1ia)
 * - F_kia = G_ia * exp(D_kia + D_1ia), for k=2..K-1
 */

class DiffSelFitnessArray : public SimpleBidimArray<vector<double>> {
  public:
    //! constructor, parameterized by baseline (G), delta (D) and Nlevel
    DiffSelFitnessArray(const Selector<vector<double>> &inbaseline,
                        const BidimSelector<vector<double>> &indelta,
                        int inNlevel)
        :  // DiffSelFitnessArray(const Selector<vector<double> >& inbaseline,
           // const BidimSelector<vector<double> >& indelta, const
           // vector<vector<int> >& inpattern) :
          SimpleBidimArray<vector<double>>(indelta.GetNrow() + 1, indelta.GetNcol(),
                                           vector<double>(indelta.GetVal(0, 0).size(), 0)),
          baseline(inbaseline),
          delta(indelta),
          Nlevel(inNlevel) {
        // baseline(inbaseline), delta(indelta), pattern(inpattern) {
        Update();
    }

    //! returns dimension of fitness profiles (should normally be 20)
    int GetDim() const { return GetVal(0, 0).size(); }

    //! full update of the array
    void Update() {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                Update(i, j);
            }
        }
    }

    //! update of column j (i.e. site j)
    void UpdateColumn(int j) {
        for (int i = 0; i < GetNrow(); i++) {
            Update(i, j);
        }
    }

    //! update of column j (i.e. site j), but only for those conditions that are
    //! flagged
    void UpdateColumn(int j, const vector<int> &flag) {
        for (int i = 0; i < GetNrow(); i++) {
            if (flag[i]) {
                Update(i, j);
            }
        }
    }

    //! update of column j (i.e. site j) and condition i
    void Update(int i, int j) {
        vector<double> &x = (*this)(i, j);
        double total = 0;
        for (int k = 0; k < GetDim(); k++) {
            double d = 0;
            if ((i == 1) || (Nlevel == 2)) {
                d += delta.GetVal(0, j)[k];
            }
            if (i > 1) {
                d += delta.GetVal(i - 1, j)[k];
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
        for (int k = 0; k < GetDim(); k++) {
            x[k] /= total;
        }
    }

  protected:
    const Selector<vector<double>> &baseline;
    const BidimSelector<vector<double>> &delta;
    // const vector<vector<int> >& pattern;
    int Nlevel;
};

#endif
