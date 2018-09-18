#include <cmath>
#include <fstream>
#include "DiffSelDoublySparseModel.hpp"
#include "Sample.hpp"
using namespace std;

/**
 * \brief An MCMC sample for DiffSelDoublySparseModel
 *
 * implements a read function, returning the MCMC estimate of the posterior
 * probabilities of differential effects across sites
 */

class DiffSelDoublySparseSample : public Sample {
  private:
    string modeltype;
    string datafile;
    string treefile;
    int ncond, nlevel, codonmodel;
    double fitnessshape;
    int fitnesscentermode;
    double epsilon;
    double pihypermean;
    double shiftprobmean;
    double shiftprobinvconc;
    int chainburnin;

  public:
    string GetModelType() override { return modeltype; }

    DiffSelDoublySparseModel *GetModel() { return (DiffSelDoublySparseModel *)model; }

    //! \brief Constructor (file name, burn-in, thinning and upper limit, see
    //! Sample)
    DiffSelDoublySparseSample(string filename, int inburnin, int inevery, int inuntil)
        : Sample(filename, inburnin, inevery, inuntil) {
        Open();
    }

    void Open() override {
        // open <name>.param
        ifstream is((name + ".param").c_str());

        // check that file exists
        if (!is) {
            cerr << "error : cannot find file : " << name << ".param\n";
            exit(1);
        }

        // read model type, and other standard fields
        is >> modeltype;
        is >> datafile >> treefile >> ncond >> nlevel;
        is >> codonmodel;
        is >> fitnessshape;
        is >> fitnesscentermode;
        is >> epsilon;
        is >> pihypermean >> shiftprobmean >> shiftprobinvconc;
        int check;
        is >> check;
        if (check) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainburnin;
        is >> chainevery >> chainuntil >> chainsaveall >> chainsize;

        if (burnin < chainburnin) {
            cerr << "error: sample burnin smaller than chain burnin\n";
            exit(1);
        }

        // make a new model depending on the type obtained from the file
        if (modeltype == "DIFFSELSPARSE") {
            model = new DiffSelDoublySparseModel(datafile, treefile, ncond, nlevel, codonmodel,
                epsilon, fitnessshape, pihypermean, shiftprobmean, shiftprobinvconc);
            GetModel()->SetFitnessCenterMode(fitnesscentermode);
        } else {
            cerr << "error when opening file " << name << '\n';
            exit(1);
        }

        GetModel()->Allocate();

        // read model (i.e. chain's last point) from <name>.param
        model->FromStream(is);

        // open <name>.chain, and prepare stream and stream iterator
        OpenChainFile();
        // now, size is defined (it is the total number of points with which this
        // Sample object will make all its various posterior averages) all these
        // points can be accessed to (only once) by repeated calls to GetNextPoint()
    }

    //! \brief computes the posterior mean estimate (and the posterior standard
    //! deviation) of omega
    void ReadPP(double cutoff, int siteoffset) {
        int Nsite = GetModel()->GetNsite();
        int Ncond = ncond;
        vector<vector<vector<double>>> pp(
            Ncond - 1, vector<vector<double>>(Nsite, vector<double>(Naa, 0)));
        vector<vector<double>> sitepp(Ncond - 1, vector<double>(Nsite, 0));
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            const vector<vector<int>> &mask = GetModel()->GetMaskArray();
            for (int k = 1; k < Ncond; k++) {
                const vector<vector<int>> &toggle = GetModel()->GetCondToggleArray(k);
                for (int j = 0; j < Nsite; j++) {
                    int n = 0;
                    for (int a = 0; a < Naa; a++) {
                        pp[k - 1][j][a] += toggle[j][a] * mask[j][a];
                        n += toggle[j][a] * mask[j][a];
                    }
                    if (n) { sitepp[k - 1][j]++; }
                }
            }
        }
        cerr << '\n';

        // normalization
        for (int k = 1; k < Ncond; k++) {
            for (int j = 0; j < Nsite; j++) {
                for (int a = 0; a < Naa; a++) { pp[k - 1][j][a] /= size; }
                sitepp[k - 1][j] /= size;
            }
        }

        // write output
        for (int k = 1; k < Ncond; k++) {
            ostringstream s;
            s << name << "_" << k << ".sitepp";
            ofstream os(s.str().c_str());
            os << "site\tsitepp";
            for (int a = 0; a < Naa; a++) { os << '\t' << AminoAcids[a]; }
            os << '\n';
            for (int j = 0; j < Nsite; j++) {
                if (sitepp[k - 1][j] > cutoff) {
                    os << j + siteoffset << '\t' << int(100 * sitepp[k - 1][j]);
                    for (int a = 0; a < Naa; a++) { os << '\t' << int(100 * pp[k - 1][j][a]); }
                    os << '\n';
                }
            }
        }

        // write output
        for (int k = 1; k < Ncond; k++) {
            ostringstream s;
            s << name << "_" << k << ".pp";
            ofstream os(s.str().c_str());
            os << "site\tAA\tpp\n";
            for (int j = 0; j < Nsite; j++) {
                for (int a = 0; a < Naa; a++) {
                    if (pp[k - 1][j][a] > cutoff) {
                        os << j + siteoffset << '\t' << AminoAcids[a] << '\t'
                           << int(100 * pp[k - 1][j][a]) << '\n';
                    }
                }
            }
        }

        cerr << "post probs written in .pp files\n";
        cerr << '\n';
    }
};

int main(int argc, char *argv[]) {
    int burnin = 0;
    int every = 1;
    int until = -1;
    string name;
    int ppred = 0;

    int siteoffset = 0;
    double cutoff = 0.90;

    try {
        if (argc == 1) { throw(0); }

        int i = 1;
        while (i < argc) {
            string s = argv[i];
            if (s == "-ppred") {
                ppred = 1;
            } else if ((s == "-x") || (s == "-extract")) {
                i++;
                if (i == argc) throw(0);
                s = argv[i];
                if (!IsInt(s)) { throw(0); }
                burnin = atoi(argv[i]);
                i++;
                if (i == argc) throw(0);
                s = argv[i];
                if (IsInt(s)) {
                    every = atoi(argv[i]);
                    i++;
                    if (i == argc) throw(0);
                    s = argv[i];
                    if (IsInt(s)) {
                        until = atoi(argv[i]);
                    } else {
                        i--;
                    }
                } else {
                    i--;
                }
            } else if (s == "-c") {
                i++;
                cutoff = atof(argv[i]);
            } else if (s == "-s") {
                i++;
                siteoffset = atoi(argv[i]);
            } else {
                if (i != (argc - 1)) { throw(0); }
                name = argv[i];
            }
            i++;
        }
        if (name == "") { throw(0); }
    } catch (...) {
        cerr << "readglobom [-x <burnin> <every> <until>] <chainname> \n";
        cerr << '\n';
        exit(1);
    }

    DiffSelDoublySparseSample *sample = new DiffSelDoublySparseSample(name, burnin, every, until);
    if (ppred) {
        sample->PostPred();
    } else {
        sample->ReadPP(cutoff, siteoffset);
    }
}
