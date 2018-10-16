#pragma once

#include <map>
#include "Array.hpp"
#include "BidimArray.hpp"
#include "BranchArray.hpp"
#include "NodeArray.hpp"
#include "PhyloProcess.hpp"
#include "PoissonRandomField.hpp"
#include "SubMatrix.hpp"
#include "SuffStat.hpp"

/**
 * \brief A general sufficient statistic for substitution histories, as a
 * function of the substitution rate matrix.
 *
 * The probability of a set of detailed substitution histories (collectively
 * denoted as S), as a function of some rate matrix Q = (Q_ab), with equilibrium
 * frequencies pi = (pi_a), can be written as:
 *
 * p(S | Q) propto (prod_a pi_a^u_a) (prod_a exp(t_a Q_aa)) (prod_ab Q_ab^v_ab),
 *
 * where u_a is the total number of times state a was seen at the root (root
 * count statistic), v_ab (pair is the total number of substitution events from
 * a to b (pair count stat), and t_a is the total waiting time in state a
 * (waiting time stat) -- all this, across all substitution histories included
 * in S.
 *
 * PolySuffStat implements this idea, by providing methods for gathering
 * sufficient statistics across substitution histories (see also
 * BranchSitePoly::AddPolySuffStat), adding them across sites and/or branches,
 * and calculating the log p(S | Q) for any matrix Q
 *
 * These Poly suffstats can be used for any Markovian substitution process (any
 * Q). In some cases (i.e. for Muse and Gaut codon models), they can be
 * furthered simplified, as a function of the nucleotide rate parameters or the
 * omega parameter of the Q matrix, leading to even more compact suff stats (see
 * OmegaPolySuffStat and NucPolySuffStat).
 *
 * In terms of implementation, these suffstats are encoded as sparse data
 * structures (since a very small subset of all possible pairs of codons will
 * typcially be visited by the substitution history of a given site, for
 * instance). This sparse encoding is crucial for efficiency (both in terms of
 * time and in terms of RAM usage).
 */

class PolySuffStat : public SuffStat {
  public:
    PolySuffStat() {}

    ~PolySuffStat() {}

    //! set suff stats to 0
    void Clear() { polycount.clear(); }

    void IncrementPolyCount(std::tuple<int, int, unsigned, unsigned> poly_tuple) {
        polycount[poly_tuple]++;
    }

    void AddPairCount(std::tuple<int, int, unsigned, unsigned> poly_tuple, int in) {
        polycount[poly_tuple] += in;
    }

    //! add Poly sufficient statistics from PhyloProcess (site-homogeneous case)
    void AddSuffStat(const PhyloProcess &process) { process.AddPolySuffStat(*this); }

    void Add(const PolySuffStat &suffstat) {
        for (auto const &i : suffstat.GetPairCountMap()) { AddPairCount(i.first, i.second); }
    }

    PolySuffStat &operator+=(const PolySuffStat &from) {
        Add(from);
        return *this;
    }

    int GetPairCount(std::tuple<int, int, unsigned, unsigned> poly_tuple) const {
        auto const i = polycount.find(poly_tuple);
        if (i == polycount.end()) { return 0; }
        return i->second;
    }

    //! return log p(S | Q) as a function of the fitness vector, the nucmatrix and theta
    double GetLogProb(PoissonRandomField &poissonrandomfield,
        const std::vector<double> &aafitnessarray, const GTRSubMatrix &nucmatrix,
        const double &theta) const {
        double total = 0;
        for (auto const &i : polycount) {
            double proba = poissonrandomfield.GetProb(std::get<0>(i.first), std::get<1>(i.first),
                std::get<2>(i.first), std::get<3>(i.first), &aafitnessarray, &nucmatrix, &theta);
            if (proba > 0) {
                total += i.second * log(proba);
            } else {
                total = -std::numeric_limits<double>::infinity();
            }
        }
        return total;
    }

    //! const access to the ordered map giving the pair count stat (sparse data
    //! structure)
    const std::map<std::tuple<int, int, unsigned, unsigned>, int> &GetPairCountMap() const {
        return polycount;
    }

  private:
    std::map<std::tuple<int, int, unsigned, unsigned>, int> polycount;
};

/**
 * \brief An array of substitution Poly sufficient statistics
 *
 * This class provides an interface for dealing with cases where
 * each item (each site, or each component of a mixture) has a different rate
 * matrix Q_i
 */

class PolySuffStatArray : public SimpleArray<PolySuffStat> {
  public:
    PolySuffStatArray(int insize) : SimpleArray<PolySuffStat>(insize) {}

    ~PolySuffStatArray() {}

    //! set all suff stats to 0
    void Clear() {
        for (int i = 0; i < GetSize(); i++) { (*this)[i].Clear(); }
    }

    //! add Poly sufficient statistics from PhyloProcess (site-heterogeneous case)
    void AddSuffStat(const PhyloProcess &process) { process.AddPolySuffStat(*this); }

    //! return total log prob (summed over all items), given an array of rate
    //! matrices
    double GetLogProb(PoissonRandomField &poissonrandomfield,
        const Selector<std::vector<double>> &siteaafitnessarray, const GTRSubMatrix &nucmatrix,
        const double &theta) const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) {
            total += GetVal(i).GetLogProb(
                poissonrandomfield, siteaafitnessarray.GetVal(i), nucmatrix, theta);
        }
        return total;
    }

    //! \brief add suffstatarray given as argument to this array based on the
    //! allocations provided as the second argument (mixture models)
    //!
    //! specifically, for each i=0..GetSize()-1, (*this)[alloc[i]] +=
    //! suffstatarray[i]
    void Add(const Selector<PolySuffStat> &suffstatarray, const Selector<int> &alloc) {
        for (int i = 0; i < suffstatarray.GetSize(); i++) {
            (*this)[alloc.GetVal(i)] += suffstatarray.GetVal(i);
        }
    }
};
