#include <cmath>
#include <fstream>
#include "AAMutSelSparseOmegaModel.hpp"
#include "Sample.hpp"
using namespace std;

/**
 * \brief An MCMC sample for AAMutSelSparseOmegaModel
 *
 * implements a read function, returning the MCMC estimate of the posterior
 * probabilities of differential effects across sites
 */

class AAMutSelSparseOmegaSample : public Sample {
  private:
    string modeltype;
    string datafile;
    string treefile;
    int omegamode, omegaprior;
    double dposompi, dposomhypermean, dposomhyperinvshape;
    double maxdposom;
    int fixhyper;
    int fixfitness;
    // -1: estimated
    double epsilon;
    double pi;

  public:
    string GetModelType() override { return modeltype; }

    AAMutSelSparseOmegaModel *GetModel() override { return (AAMutSelSparseOmegaModel *)model; }

    //! \brief Constructor (file name, burn-in, thinning and upper limit, see
    //! Sample)
    AAMutSelSparseOmegaSample(string filename, int inburnin, int inevery, int inuntil)
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
        is >> datafile >> treefile;
        is >> omegamode >> omegaprior >> dposompi >> dposomhypermean >> dposomhyperinvshape >> maxdposom;
        is >> fixhyper;
        is >> fixfitness;
        is >> epsilon;
        is >> pi;
        int check;
        is >> check;
        if (check) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "AAMUTSELSparseOMEGA") {
            model = new AAMutSelSparseOmegaModel(datafile,treefile,omegamode,omegaprior,fixhyper,fixfitness,epsilon,pi);

            if (omegaprior != 0) {
                GetModel()->SetMaxDPosOm(maxdposom);
                GetModel()->SetDPosOmHyperParameters(dposompi, dposomhypermean,
                                                     dposomhyperinvshape);
            }
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
};

int main(int argc, char *argv[]) {
    int burnin = 0;
    int every = 1;
    int until = -1;
    string name;
    int ppred = 0;

    try {
        if (argc == 1) {
            throw(0);
        }

        int i = 1;
        while (i < argc) {
            string s = argv[i];
            if (s == "-ppred") {
                ppred = 1;
            } else if (s == "-allppred")    {
                ppred = 2;
            } else if ((s == "-x") || (s == "-extract")) {
                i++;
                if (i == argc) throw(0);
                s = argv[i];
                if (!IsInt(s)) {
                    throw(0);
                }
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
            } else {
                if (i != (argc - 1)) {
                    throw(0);
                }
                name = argv[i];
            }
            i++;
        }
        if (name == "") {
            throw(0);
        }
    } catch (...) {
        cerr << "readglobom [-x <burnin> <every> <until>] <chainname> \n";
        cerr << '\n';
        exit(1);
    }

    AAMutSelSparseOmegaSample *sample = new AAMutSelSparseOmegaSample(name, burnin, every, until);
    if (ppred == 2) {
        sample->AllPostPred();
    } else if (ppred == 1)  {
        sample->PostPred();
    } else  {
        cerr << "read not yet implemented\n";
        exit(1);
    }
}
