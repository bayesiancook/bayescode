
#ifndef GTRSUFFSTAT_H
#define GTRSUFFSTAT_H

#include "AASubSelSubMatrix.hpp"
#include "AASubSelSubMatrixArray.hpp"
#include "GTRSubMatrix.hpp"
#include "GTRSubMatrixArray.hpp"
#include "PathSuffStat.hpp"

class RelRateSuffStat : public SuffStat {
  public:
    RelRateSuffStat(int innstate) : nstate(innstate) {
        nrr = innstate * (innstate - 1) / 2;
        rrcount.assign(nrr, 0);
        rrbeta.assign(nrr, 0);
    }

    ~RelRateSuffStat() {}

    int GetCount(int i) const { return rrcount[i]; }

    double GetBeta(int i) const { return rrbeta[i]; }

    void Clear() {
        for (int i = 0; i < nrr; i++) {
            rrcount[i] = 0;
            rrbeta[i] = 0;
        }
    }

    void AddSuffStat(
        const GTRSubMatrix &matrix, const PathSuffStat &pathsuffstat, double rate = 1) {
        const std::map<int, double> &waitingtime = pathsuffstat.GetWaitingTimeMap();
        for (std::map<int, double>::const_iterator i = waitingtime.begin(); i != waitingtime.end();
             i++) {
            int a = i->first;
            for (int b = 0; b < nstate; b++) {
                if (b != a) { rrbeta[rrindex(a, b)] += matrix.Stationary(b) * i->second * rate; }
            }
        }

        const std::map<pair<int, int>, int> &paircount = pathsuffstat.GetPairCountMap();
        for (std::map<pair<int, int>, int>::const_iterator i = paircount.begin();
             i != paircount.end(); i++) {
            rrcount[rrindex(i->first.first, i->first.second)] += i->second;
        }
    }

    void AddSuffStat(
        const GTRSubMatrixArray &matrixarray, const PathSuffStatArray &pathsuffstatarray) {
        for (int i = 0; i < matrixarray.GetSize(); i++) {
            AddSuffStat(matrixarray.GetVal(i), pathsuffstatarray.GetVal(i));
        }
    }

    void AddSuffStat(const GTRSubMatrixArray &matrixarray,
        const PathSuffStatArray &pathsuffstatarray, const Selector<double> &rate) {
        for (int i = 0; i < matrixarray.GetSize(); i++) {
            AddSuffStat(matrixarray.GetVal(i), pathsuffstatarray.GetVal(i), rate.GetVal(i));
        }
    }

    void AddSuffStat(
        const AASubSelSubMatrix &matrix, const PathSuffStat &pathsuffstat, double rate = 1) {
        const std::map<int, double> &waitingtime = pathsuffstat.GetWaitingTimeMap();
        for (std::map<int, double>::const_iterator i = waitingtime.begin(); i != waitingtime.end();
             i++) {
            int a = i->first;
            for (int b = 0; b < nstate; b++) {
                if (b != a) {
                    rrbeta[rrindex(a, b)] +=
                        matrix(a, b) / matrix.RelativeRate(a, b) * i->second * rate;
                }
            }
        }

        const std::map<pair<int, int>, int> &paircount = pathsuffstat.GetPairCountMap();
        for (std::map<pair<int, int>, int>::const_iterator i = paircount.begin();
             i != paircount.end(); i++) {
            rrcount[rrindex(i->first.first, i->first.second)] += i->second;
        }
    }

    void AddSuffStat(
        const AASubSelSubMatrixArray &matrixarray, const PathSuffStatArray &pathsuffstatarray) {
        for (int i = 0; i < matrixarray.GetSize(); i++) {
            AddSuffStat(matrixarray.GetVal(i), pathsuffstatarray.GetVal(i));
        }
    }

    void AddSuffStat(const AASubSelSubMatrixArray &matrixarray,
        const PathSuffStatArray &pathsuffstatarray, const Selector<double> &rate) {
        for (int i = 0; i < matrixarray.GetSize(); i++) {
            AddSuffStat(matrixarray.GetVal(i), pathsuffstatarray.GetVal(i), rate.GetVal(i));
        }
    }

    double GetLogProb(const vector<double> &rr) const {
        double total = 0;
        for (int i = 0; i < nrr; i++) { total += rrcount[i] * log(rr[i]) - rrbeta[i] * rr[i]; }
        return total;
    }

    void Add(const RelRateSuffStat &from) {
        for (int i = 0; i < nrr; i++) {
            rrcount[i] += from.rrcount[i];
            rrbeta[i] += from.rrbeta[i];
        }
    }

    RelRateSuffStat &operator+=(const RelRateSuffStat &from) {
        Add(from);
        return *this;
    }

  private:
    int rrindex(int i, int j) {
        return (i < j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1
                       : (2 * nstate - j - 1) * j / 2 + i - j - 1;
    }

    int nstate;
    int nrr;
    vector<int> rrcount;
    vector<double> rrbeta;
};

class ProfileSuffStat : public SuffStat {
  public:
    ProfileSuffStat(int innstate)
        : nstate(innstate), profilecount(innstate, 0), profilebeta(innstate, 0) {}
    ~ProfileSuffStat() {}

    void Clear() {
        for (int i = 0; i < nstate; i++) {
            profilecount[i] = 0;
            profilebeta[i] = 0;
        }
    }

    void AddSuffStat(
        const GTRSubMatrix &matrix, const PathSuffStat &pathsuffstat, double rate = 1) {
        const std::map<int, int> &rootcount = pathsuffstat.GetRootCountMap();
        for (std::map<int, int>::const_iterator i = rootcount.begin(); i != rootcount.end(); i++) {
            profilecount[i->first] += i->second;
        }

        const std::map<int, double> &waitingtime = pathsuffstat.GetWaitingTimeMap();
        for (std::map<int, double>::const_iterator i = waitingtime.begin(); i != waitingtime.end();
             i++) {
            int a = i->first;
            for (int b = 0; b < nstate; b++) {
                if (b != a) { profilebeta[b] += matrix.RelativeRate(a, b) * i->second * rate; }
            }
        }

        const std::map<pair<int, int>, int> &paircount = pathsuffstat.GetPairCountMap();
        for (std::map<pair<int, int>, int>::const_iterator i = paircount.begin();
             i != paircount.end(); i++) {
            profilecount[i->first.second] += i->second;
        }
    }

    void AddSuffStat(
        const GTRSubMatrixArray &matrixarray, const PathSuffStatArray &pathsuffstatarray) {
        for (int i = 0; i < matrixarray.GetSize(); i++) {
            AddSuffStat(matrixarray.GetVal(i), pathsuffstatarray.GetVal(i));
        }
    }

    void AddSuffStat(const GTRSubMatrixArray &matrixarray,
        const PathSuffStatArray &pathsuffstatarray, const Selector<double> &rate) {
        for (int i = 0; i < matrixarray.GetSize(); i++) {
            AddSuffStat(matrixarray.GetVal(i), pathsuffstatarray.GetVal(i), rate.GetVal(i));
        }
    }

    double GetLogProb(const vector<double> &profile) const {
        double total = 0;
        for (int i = 0; i < nstate; i++) {
            total += profilecount[i] * log(profile[i]) - profilebeta[i] * profile[i];
        }
        return total;
    }

    void Add(const ProfileSuffStat &from) {
        for (int i = 0; i < nstate; i++) {
            profilecount[i] += from.profilecount[i];
            profilebeta[i] += from.profilebeta[i];
        }
    }

    ProfileSuffStat &operator+=(const ProfileSuffStat &from) {
        Add(from);
        return *this;
    }

  private:
    int nstate;
    vector<int> profilecount;
    vector<double> profilebeta;
};

class ProfileSuffStatArray : public SimpleArray<ProfileSuffStat> {
  public:
    ProfileSuffStatArray(int insize, int innstate)
        : SimpleArray<ProfileSuffStat>(insize, ProfileSuffStat(innstate)){};
    ~ProfileSuffStatArray() {}

    void Clear() {
        for (int i = 0; i < GetSize(); i++) { (*this)[i].Clear(); }
    }

    void AddSuffStat(const Selector<GTRSubMatrix> &matrixarray,
        const Selector<PathSuffStat> &pathsuffstatarray) {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i].AddSuffStat(matrixarray.GetVal(i), pathsuffstatarray.GetVal(i));
        }
    }

    void AddSuffStat(const Selector<GTRSubMatrix> &matrixarray,
        const Selector<PathSuffStat> &pathsuffstatarray, const Selector<double> &rate) {
        for (int i = 0; i < GetSize(); i++) {
            (*this)[i].AddSuffStat(
                matrixarray.GetVal(i), pathsuffstatarray.GetVal(i), rate.GetVal(i));
        }
    }

    double GetLogProb(const Selector<vector<double>> &profilearray) const {
        double total = 0;
        for (int i = 0; i < GetSize(); i++) {
            total += GetVal(i).GetLogProb(profilearray.GetVal(i));
        }
        return total;
    }
};

#endif
