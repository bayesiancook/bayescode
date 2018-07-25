#ifndef DISCGAMMA_H
#define DISCGAMMA_H

#include "Array.hpp"
#include "IncompleteGamma.hpp"

class DiscGamma : public SimpleArray<double> {
  public:
    DiscGamma(int inncat, double inmean, double ininvshape)
        : SimpleArray<double>(inncat),
          ncat(inncat),
          mean(inmean),
          invshape(ininvshape),
          weight(inncat, 1.0 / inncat) {
        ComputeDiscGamma();
    }

    ~DiscGamma() {}

    int GetNcat() const { return ncat; }

    void SetParameters(double inmean, double ininvshape) {
        mean = inmean;
        invshape = ininvshape;
        ComputeDiscGamma();
    }

    const vector<double> &GetWeights() const { return weight; }

    double GetPostProbArray(OmegaPathSuffStatArray &suffstatarray,
                            vector<vector<double>> &postprobarray) const {
        double total = 0;
        for (int i = 0; i < GetNcat(); i++) {
            total += GetPostProbArray(suffstatarray.GetVal(i), postprobarray[i]);
        }
        return total;
    }

    double GetPostProbArray(const OmegaPathSuffStat &suffstat, vector<double> &postprob) const {
        double logp[GetNcat()];
        double max = 0;
        for (int cat = 0; cat < GetNcat(); cat++) {
            logp[cat] = suffstat.GetLogProb(GetVal(cat));
            if ((!cat) || (max < logp[cat])) {
                max = logp[cat];
            }
        }
        double tot = 0;
        for (int cat = 0; cat < GetNcat(); cat++) {
            postprob[cat] = weight[cat] * exp(logp[cat] - max);
            tot += postprob[cat];
        }
        for (int cat = 0; cat < GetNcat(); cat++) {
            postprob[cat] /= tot;
        }
        double ret = log(tot) + max;
        if (std::isinf(ret)) {
            cerr << "ret is inf: " << tot << '\t' << max << '\n';
            exit(1);
        }
        if (std::isnan(ret)) {
            cerr << "ret is nan: " << tot << '\t' << max << '\n';
            for (int cat = 0; cat < GetNcat(); cat++) {
                cerr << GetVal(cat) << '\t' << weight[cat] << '\t' << logp[cat] << '\t'
                     << postprob[cat] << '\n';
            }
            cerr << tot << '\t' << log(tot) << '\t' << max << '\n';
            exit(1);
        }
        return ret;
    }

  private:
    void ComputeDiscGamma() {
        double alpha = 1.0 / invshape;

        double x[GetNcat()];
        double y[GetNcat()];
        double lg = Random::logGamma(alpha + 1.0);
        for (int i = 0; i < GetNcat(); i++) {
            x[i] = PointGamma((i + 1.0) / GetNcat(), alpha, alpha);
        }
        for (int i = 0; i < GetNcat() - 1; i++) {
            y[i] = IncompleteGamma(alpha * x[i], alpha + 1, lg);
        }
        y[GetNcat() - 1] = 1.0;
        (*this)[0] = GetNcat() * y[0];
        for (int i = 1; i < GetNcat(); i++) {
            (*this)[i] = GetNcat() * (y[i] - y[i - 1]);
        }
        for (int i = 0; i < GetNcat(); i++) {
            (*this)[i] *= mean;
        }
    }

    int ncat;
    double mean;
    double invshape;
    vector<double> weight;
};

#endif
