
#pragma once

#include <map>
#include "Array.hpp"
#include "BidimArray.hpp"
#include "BranchArray.hpp"
#include "NodeArray.hpp"
#include "SubMatrix.hpp"
#include "SuffStat.hpp"
#include "PhyloProcess.hpp"
#include "PathSuffStat.hpp"

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
 * RelativePathSuffStat implements this idea, by providing methods for gathering
 * sufficient statistics across substitution histories (see also
 * BranchSitePath::AddRelativePathSuffStat), adding them across sites and/or branches,
 * and calculating the log p(S | Q) for any matrix Q
 *
 * These path suffstats can be used for any Markovian substitution process (any
 * Q). In some cases (i.e. for Muse and Gaut codon models), they can be
 * furthered simplified, as a function of the nucleotide rate parameters or the
 * omega parameter of the Q matrix, leading to even more compact suff stats (see
 * OmegaRelativePathSuffStat and NucRelativePathSuffStat).
 *
 * In terms of implementation, these suffstats are encoded as sparse data
 * structures (since a very small subset of all possible pairs of codons will
 * typcially be visited by the substitution history of a given site, for
 * instance). This sparse encoding is crucial for efficiency (both in terms of
 * time and in terms of RAM usage).
 */

class RelativePathSuffStat {
// class RelativePathSuffStat : public SuffStat {
  public:

    RelativePathSuffStat(int inNstate) : Nstate(inNstate) {}

    RelativePathSuffStat(const RelativePathSuffStat& from) : Nstate(from.Nstate)  {
        Clear();
        Add(from);
    }

    ~RelativePathSuffStat() {}

    //! set suff stats to 0
    void Clear() {
        rootcount.clear();
        paircount.clear();
        waitingtime.clear();
    }

    void Add(const PathSuffStat &suffstat, double length) {
        for (std::map<int, double>::const_iterator i = suffstat.GetRootCountMap().begin();
             i != suffstat.GetRootCountMap().end(); i++) {
            rootcount[i->first] += i->second;
        }
        for (std::map<pair<int, int>, double>::const_iterator i = suffstat.GetPairCountMap().begin();
             i != suffstat.GetPairCountMap().end(); i++) {
            paircount[pair<int,int>(i->first.first, i->first.second)] += i->second;
        }
        for (std::map<int, double>::const_iterator i = suffstat.GetWaitingTimeMap().begin();
             i != suffstat.GetWaitingTimeMap().end(); i++) {
            waitingtime[i->first] += i->second / length;
        }
    }

    void Add(const RelativePathSuffStat &suffstat)  {
        for (std::map<int, double>::const_iterator i = suffstat.GetRootCountMap().begin();
             i != suffstat.GetRootCountMap().end(); i++) {
            rootcount[i->first] += i->second;
        }
        for (std::map<pair<int, int>, double>::const_iterator i = suffstat.GetPairCountMap().begin();
             i != suffstat.GetPairCountMap().end(); i++) {
            paircount[pair<int,int>(i->first.first, i->first.second)] += i->second;
        }
        for (std::map<int, double>::const_iterator i = suffstat.GetWaitingTimeMap().begin();
             i != suffstat.GetWaitingTimeMap().end(); i++) {
            waitingtime[i->first] += i->second;
        }
    }

    RelativePathSuffStat &operator+=(const RelativePathSuffStat &from) {
        Add(from);
        return *this;
    }

    double GetRootCount(int state) const {
        std::map<int, double >::const_iterator i = rootcount.find(state);
        if (i == rootcount.end()) {
            return 0;
        }
        return i->second;
    }

    double GetPairCount(int state1, int state2) const {
        std::map<pair<int, int>, double >::const_iterator i =
            paircount.find(pair<int, int>(state1, state2));
        if (i == paircount.end()) {
            return 0;
        }
        return i->second;
    }

    double GetWaitingTime(int state) const {
        std::map<int, double>::const_iterator i = waitingtime.find(state);
        if (i == waitingtime.end()) {
            return 0;
        }
        return i->second;
    }

    //! return log p(S | Q) as a function of the Q matrix given as the argument
    double GetLogProb(const SubMatrix &mat, double length) const {
        double total = 0;
        auto stat = mat.GetStationary();
        for (std::map<int, double >::const_iterator i = rootcount.begin(); i != rootcount.end(); i++) {
            total += i->second * log(stat[i->first]);
        }
        for (std::map<int, double>::const_iterator i = waitingtime.begin(); i != waitingtime.end();
             i++) {
            total += length * i->second * mat(i->first, i->first);
        }
        for (std::map<pair<int, int>, double >::const_iterator i = paircount.begin();
             i != paircount.end(); i++) {
            total += i->second * log(mat(i->first.first, i->first.second));
        }
        return total;
    }

    //! const access to the ordered map giving the root count stat (sparse data
    //! structure)
    const std::map<int, double > &GetRootCountMap() const { return rootcount; }
    //! const access to the ordered map giving the pair count stat (sparse data
    //! structure)
    const std::map<pair<int, int>, double > &GetPairCountMap() const { return paircount; }
    //! const access to the ordered map giving the waiting time stat (sparse data
    //! structure)
    const std::map<int, double> &GetWaitingTimeMap() const { return waitingtime; }

    //! return size of object, when put into an MPI buffer
    unsigned int GetMPISize() const {
        if (! Nstate)   {
            cerr << "error in RelativePathSuffStat::GetMPISize: Nstate not initialized\n";
            exit(1);
        }
        return 2*Nstate + Nstate*(Nstate-1);
    }

    //! put object into MPI buffer
    void MPIPut(MPIBuffer &buffer) const {
        for (int i=0; i<Nstate; i++)    {
            buffer << GetRootCount(i);
        }
        for (int i=0; i<Nstate; i++)    {
            for (int j=0; j<Nstate; j++)    {
                if (i != j) {
                    buffer << GetPairCount(i,j);
                }
            }
        }
        for (int i=0; i<Nstate; i++)    {
            buffer << GetWaitingTime(i);
        }
    }

    //! get object from MPI buffer
    void MPIGet(const MPIBuffer &buffer) {
        Clear();
        Add(buffer);
    }

    //! get a nucpath suffstat from MPI buffer and add it to this
    void Add(const MPIBuffer &buffer) {
        for (int i=0; i<Nstate; i++)    {
            double tmp;
            buffer >> tmp;
            if (tmp)    {
                rootcount[i] += tmp;
            }
        }
        for (int i=0; i<Nstate; i++)    {
            for (int j=0; j<Nstate; j++)    {
                if (i != j) {
                    double tmp;
                    buffer >> tmp;
                    if (tmp)    {
                        paircount[pair<int,int>(i,j)] += tmp;
                    }
                }
            }
        }
        for (int i=0; i<Nstate; i++)    {
            double tmp;
            buffer >> tmp;
            if (tmp)    {
                waitingtime[i] += tmp;
            }
        }
    }

    void ToStream(ostream& os) const { 
        for (int i=0; i<Nstate; i++)    {
            os << GetRootCount(i) << '\t';
        }
        for (int i=0; i<Nstate; i++)    {
            for (int j=0; j<Nstate; j++)    {
                if (i != j) {
                    os << GetPairCount(i,j) << '\t';
                }
            }
        }
        for (int i=0; i<Nstate; i++)    {
            os << GetWaitingTime(i);
            if (i < Nstate-1)   {
                os << '\t';
            }
        }
    }

    void FromStream(istream& is) { 
        Clear();
        for (int i=0; i<Nstate; i++)    {
            double tmp;
            is >> tmp;
            if (tmp)    {
                rootcount[i] = tmp;
            }
        }
        for (int i=0; i<Nstate; i++)    {
            for (int j=0; j<Nstate; j++)    {
                if (i != j) {
                    double tmp;
                    is >> tmp;
                    if (tmp)    {
                        paircount[pair<int,int>(i,j)] = tmp;
                    }
                }
            }
        }
        for (int i=0; i<Nstate; i++)    {
            double tmp;
            is >> tmp;
            if (tmp)    {
                waitingtime[i] = tmp;
            }
        }
    }

  private:
    int Nstate;
    std::map<int, double> rootcount;
    std::map<pair<int, int>, double> paircount;
    std::map<int, double> waitingtime;
};


ostream& operator<<(ostream& os, const RelativePathSuffStat& suffstat)    {
    suffstat.ToStream(os);
    return os;
}

istream& operator>>(istream& is, RelativePathSuffStat& suffstat)    {
    suffstat.FromStream(is);
    return is;
}
/**
 * \brief An array of substitution path sufficient statistics
 *
 * This class provides an interface for dealing with cases where
 * each item (each site, or each component of a mixture) has a different rate
 * matrix Q_i
 */

class RelativePathSuffStatArray : public SimpleArray<RelativePathSuffStat> {
  public:
    RelativePathSuffStatArray(int insize, int Nstate) : SimpleArray<RelativePathSuffStat>(insize, RelativePathSuffStat(Nstate)) {}
    ~RelativePathSuffStatArray() {}

    //! set all suff stats to 0
    void Clear() {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i].Clear();
        }
    }

    //! return total log prob (summed over all items), given an array of rate
    //! matrices
    double GetLogProb(const Selector<SubMatrix> &matrixarray, const Selector<double>& lengtharray) const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) {
            total += GetVal(i).GetLogProb(matrixarray.GetVal(i), lengtharray.GetVal(i));
        }
        return total;
    }

    //! get array from MPI buffer and add it to this array (member-wise addition)
    void Add(const MPIBuffer &buffer) {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i] += buffer;
        }
    }
};

/**
 * \brief A NodeArray of substitution path sufficient statistics
 *
 * This class provides an interface for dealing with cases where
 * each branch has a different rate matrix Q_j
 */

class RelativePathSuffStatNodeArray : public SimpleNodeArray<RelativePathSuffStat> {
  public:
    RelativePathSuffStatNodeArray(const Tree &intree, int Nstate) : SimpleNodeArray<RelativePathSuffStat>(intree, RelativePathSuffStat(Nstate)) {}
    ~RelativePathSuffStatNodeArray() {}

    //! set all suff stats to 0
    void Clear() {
        for (int i = 0; i < GetNnode(); i++) {
            (*this)[i].Clear();
        }
    }

    //! add path sufficient statistics from PhyloProcess (site-heterogeneous case)
    void AddSuffStat(const NodeSelector<PathSuffStat>& pathsuffstatarray, const BranchSelector<double>& lengtharray)    {
        RecursiveAddSuffStat(GetTree().GetRoot(), pathsuffstatarray, lengtharray);
    }

    void RecursiveAddSuffStat(const Link* from, const NodeSelector<PathSuffStat>& pathsuffstatarray, const BranchSelector<double>& lengtharray) {

        if (from->isRoot()) {
            (*this)[from->GetNode()->GetIndex()].Add(pathsuffstatarray.GetVal(from->GetNode()->GetIndex()), 0);
        }
        else    {
            (*this)[from->GetNode()->GetIndex()].Add(pathsuffstatarray.GetVal(from->GetNode()->GetIndex()), lengtharray.GetVal(from->GetBranch()->GetIndex()));
        }
        for (const Link *link=from->Next(); link!=from; link=link->Next()) {
            RecursiveAddSuffStat(link->Out(), pathsuffstatarray, lengtharray);
        }
    }

    //! return total log prob (summed over all items), given an array of rate
    //! matrices
    double GetLogProb(const BranchSelector<SubMatrix> &matrixarray,
                      const SubMatrix &rootmatrix, const BranchSelector<double>& lengtharray) const {
        double ret = RecursiveGetLogProb(GetTree().GetRoot(), matrixarray, rootmatrix, lengtharray);
        return ret;
    }

    double RecursiveGetLogProb(const Link *from, const BranchSelector<SubMatrix> &matrixarray,
                               const SubMatrix &rootmatrix, const BranchSelector<double>& lengtharray) const {
        double total = 0;
        if (from->isRoot()) {
            total += GetVal(from->GetNode()->GetIndex()).GetLogProb(rootmatrix, 0);
        } else {
            total += GetVal(from->GetNode()->GetIndex())
                         .GetLogProb(matrixarray.GetVal(from->GetBranch()->GetIndex()), lengtharray.GetVal(from->GetBranch()->GetIndex()));
        }
        for (const Link *link=from->Next(); link!=from; link=link->Next()) {
            total += RecursiveGetLogProb(link->Out(), matrixarray, rootmatrix, lengtharray);
        }
        return total;
    }

    //! get array from MPI buffer and add it to this array (member-wise addition)
    void Add(const MPIBuffer &buffer) {
        for (int i = 0; i < GetNnode(); i++) {
            (*this)[i] += buffer;
        }
    }
};

