#ifndef MECHM2AMIX_H
#define MECHM2AMIX_H

#include "Array.hpp"

class MechM2aMix : public SimpleArray<double> {
  public:
    MechM2aMix(double posom, double posw) : SimpleArray<double>(2), weight(2) {
        (*this)[0] = 1.0;
        (*this)[1] = posom;
        weight[0] = posw;
        weight[1] = (1 - posw);
    }

    ~MechM2aMix() {}

    void SetParameters(double posom, double posw) {
        (*this)[0] = 1.0;
        (*this)[1] = posom;
        weight[0] = posw;
        weight[1] = (1 - posw);
    }

    const vector<double> &GetWeights() const { return weight; }

    double GetPostProbArray(const OmegaPathSuffStat &suffstat, vector<double> &postprob) const {
        double logp[GetSize()];
        double max = 0;
        for (int cat = 0; cat < GetSize(); cat++) {
            logp[cat] = suffstat.GetLogProb(GetVal(cat));
            if ((!cat) || (max < logp[cat])) { max = logp[cat]; }
        }
        double tot = 0;
        for (int cat = 0; cat < GetSize(); cat++) {
            postprob[cat] = weight[cat] * exp(logp[cat] - max);
            tot += postprob[cat];
        }
        for (int cat = 0; cat < GetSize(); cat++) { postprob[cat] /= tot; }
        double ret = log(tot) + max;
        if (std::isinf(ret)) {
            cerr << "ret is inf: " << tot << '\t' << max << '\n';
            exit(1);
        }
        if (std::isnan(ret)) {
            cerr << "ret is nan: " << tot << '\t' << max << '\n';
            for (int cat = 0; cat < GetSize(); cat++) {
                cerr << GetVal(cat) << '\t' << weight[cat] << '\t' << logp[cat] << '\t'
                     << postprob[cat] << '\n';
            }
            cerr << tot << '\t' << log(tot) << '\t' << max << '\n';
            exit(1);
        }
        return ret;
    }

    double GetPostProbArray(
        const OmegaPathSuffStatArray &suffstatarray, vector<vector<double>> &postprobarray) const {
        double total = 0;
        for (int i = 0; i < suffstatarray.GetSize(); i++) {
            total += GetPostProbArray(suffstatarray.GetVal(i), postprobarray[i]);
        }
        return total;
    }

  private:
    vector<double> weight;
};

#endif
