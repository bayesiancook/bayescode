#pragma once

#include <cmath>
#include "Array.hpp"
#include "CodonSuffStat.hpp"

class M2aMix : public SimpleArray<double> {
  public:
    M2aMix(double purom, double posom, double purw, double posw)
        : SimpleArray<double>(3), weight(3) {
        (*this)[0] = purom;
        (*this)[1] = 1.0;
        (*this)[2] = posom;
        weight[0] = purw * (1 - posw);
        weight[1] = (1 - purw) * (1 - posw);
        weight[2] = posw;
    }

    ~M2aMix() {}

    void SetParameters(double purom, double posom, double purw, double posw) {
        (*this)[0] = purom;
        (*this)[1] = 1.0;
        (*this)[2] = posom;
        weight[0] = purw * (1 - posw);
        weight[1] = (1 - purw) * (1 - posw);
        weight[2] = posw;
    }

    const std::vector<double> &GetWeights() const { return weight; }

    double GetPostProbArray(
        const OmegaPathSuffStat &suffstat, std::vector<double> &postprob) const {
        double logp[GetSize()];
        double max = 0;
        int first = 1;
        for (int cat = 0; cat < GetSize(); cat++) {
            logp[cat] = suffstat.GetLogProb(GetVal(cat));
            if (first || (max < logp[cat])) {
                if (weight[cat]) {
                    max = logp[cat];
                    first = 0;
                }
            }
        }
        double tot = 0;
        for (int cat = 0; cat < GetSize(); cat++) {
            if (weight[cat]) {
                postprob[cat] = weight[cat] * exp(logp[cat] - max);
            } else {
                postprob[cat] = 0;
            }
            tot += postprob[cat];
        }
        for (int cat = 0; cat < GetSize(); cat++) { postprob[cat] /= tot; }
        double ret = log(tot) + max;
        if (std::isinf(ret)) {
            std::cerr << "in M2aMix::GetPostProbArray\n";
            std::cerr << "ret is inf: " << tot << '\t' << max << '\n';
            std::cerr << "omega suff stat: " << suffstat.GetCount() << '\t' << suffstat.GetBeta()
                      << '\n';
            for (int cat = 0; cat < GetSize(); cat++) {
                std::cerr << GetVal(cat) << '\t' << weight[cat] << '\t' << logp[cat] << '\t'
                          << postprob[cat] << '\n';
            }
            std::cerr << tot << '\t' << log(tot) << '\t' << max << '\n';
            exit(1);
        }
        if (std::isnan(ret)) {
            std::cerr << "ret is nan: " << tot << '\t' << max << '\n';
            for (int cat = 0; cat < GetSize(); cat++) {
                std::cerr << GetVal(cat) << '\t' << weight[cat] << '\t' << logp[cat] << '\t'
                          << postprob[cat] << '\n';
            }
            std::cerr << tot << '\t' << log(tot) << '\t' << max << '\n';
            exit(1);
        }
        return ret;
    }

    double GetPostProbArray(const OmegaPathSuffStatArray &suffstatarray,
        std::vector<std::vector<double>> &postprobarray) const {
        double total = 0;
        for (int i = 0; i < suffstatarray.GetSize(); i++) {
            total += GetPostProbArray(suffstatarray.GetVal(i), postprobarray[i]);
        }
        return total;
    }

  private:
    std::vector<double> weight;
};
