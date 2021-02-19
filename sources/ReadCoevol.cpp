
#include <cmath>
#include <fstream>
#include "Sample.hpp"
#include "CoevolModel.hpp"
#include "DistBranchArray.hpp"
using namespace std;

/**
 * \brief An MCMC sample for CoevolModel
 *
 * implements a simple read function, returning the MCMC estimate of the
 * posterior mean and standard deviation of omega=dN/dS
 */

class CoevolSample : public Sample {
  private:
    string modeltype;
    string datafile, contdatafile, treefile, rootfile;
    GeneticCodeType codetype;

  public:
    string GetModelType() override { return modeltype; }

    CoevolModel *GetModel() override { return (CoevolModel *)model; }

    //! \brief Constructor (file name, burn-in, thinning and upper limit, see
    //! Sample)
    CoevolSample(string filename, int inburnin, int inevery, int inuntil)
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
        is >> datafile >> contdatafile >> treefile >> rootfile;
        is >> codetype;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "COEVOLDNDS") {
            model = new CoevolModel(datafile, contdatafile, treefile, rootfile, codetype);
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
    void Read() {
        cerr << size << " points to read\n";

        double meanomega = 0;
        double varomega = 0;

        DistBranchArray<double> omegatree(GetModel()->GetTree());

        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->FastUpdate();
            omegatree.Add(GetModel()->GetOmegaTree());
            double om = GetModel()->GetMeanOmega();
            meanomega += om;
            varomega += om * om;
        }
        cerr << '\n';
        meanomega /= size;
        varomega /= size;
        varomega -= meanomega * meanomega;

        cout << "posterior mean omega : " << meanomega << '\t' << sqrt(varomega) << '\n';
        ofstream os((name + ".postmeanomega.tre").c_str());
        omegatree.MeanToStream(os);
    }
};

int main(int argc, char *argv[]) {
    int burnin = 0;
    int every = 1;
    int until = -1;
    int ppred = 0;

    string name;

    try {
        if (argc == 1) {
            throw(0);
        }

        int i = 1;
        while (i < argc) {
            string s = argv[i];
            if ((s == "-x") || (s == "-extract")) {
                i++;
                if (i == argc) throw(0);
                s = argv[i];
                burnin = atoi(argv[i]);
                i++;
                if (i == argc) throw(0);
                s = argv[i];
                every = atoi(argv[i]);
                i++;
                if (i == argc) throw(0);
                s = argv[i];
                until = atoi(argv[i]);
            } else if (s == "-ppred") {
                ppred = 1;
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

    CoevolSample *sample = new CoevolSample(name, burnin, every, until);
    if (ppred) {
        sample->PostPred();
    } else {
        sample->Read();
    }
}
