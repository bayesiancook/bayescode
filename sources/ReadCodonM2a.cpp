#include <cmath>
#include <fstream>
#include "CodonM2aModel.hpp"
#include "Sample.hpp"
using namespace std;

/**
 * \brief An MCMC sample for CodonM2aModel
 *
 * implements a simple read function, returning the MCMC estimate of the
 * posterior mean and standard deviation of omega=dN/dS
 */

class CodonM2aSample : public Sample {
  private:
    string modeltype;
    string datapath;
    string newpath;
    string datafile;
    string treefile;
    double pi;
    double puromhypermean, puromhyperinvconc, dposomhypermean, dposomhyperinvshape;
    double purwhypermean, purwhyperinvconc, poswhypermean, poswhyperinvconc;

  public:
    string GetModelType() override { return modeltype; }

    CodonM2aModel *GetModel() override { return (CodonM2aModel *)model; }

    //! \brief Constructor (file name, burn-in, thinning and upper limit, see
    //! Sample)
    CodonM2aSample(string innewpath, string filename, int inburnin, int inevery, int inuntil)
        : Sample(filename, inburnin, inevery, inuntil) {
        newpath = innewpath;
        Open();
    }

    void PostPredSimu(double shrinkposw, double shrinkdposom) {
        cerr << size << " points to read\n";
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->ShrinkMixtureParameters(shrinkposw, shrinkdposom);
           //  GetModel()->GetPurOm(), 0, GetModel()->GetPurW(), 0);
            ostringstream s;
            s << "ppred" << name << "_" << i;
            // s << "ppred" << name << "_" << i << ".ali";
            model->PostPred(s.str());
        }
        cerr << '\n';
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
        is >> datapath >> datafile >> treefile;
        if (newpath != "None")  {
            datapath = newpath;
        }
        is >> pi;
        is >> puromhypermean >> puromhyperinvconc;
        is >> dposomhypermean >> dposomhyperinvshape;
        is >> purwhypermean >> purwhyperinvconc;
        is >> poswhypermean >> poswhyperinvconc;
        int check;
        is >> check;
        if (check) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "CODONM2A") {
            model = new CodonM2aModel(datapath, datafile, treefile, pi);
        } else {
            cerr << "error when opening file " << name << '\n';
            exit(1);
        }

        GetModel()->SetMixtureHyperParameters(puromhypermean, puromhyperinvconc, dposomhypermean,
                                              dposomhyperinvshape, pi, purwhypermean, purwhyperinvconc,
                                              poswhypermean, poswhyperinvconc);
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
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
        }
        cerr << '\n';
    }
};

int main(int argc, char *argv[]) {
    string newpath = "None";
    int burnin = 0;
    int every = 1;
    int until = -1;
    int ppred = 0;
    double shrinkposw = 1.0;
    double shrinkdposom = 1.0;

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
            } else if (s == "-p") {
                i++;
                newpath = argv[i];
            } else if (s == "-ppred") {
                ppred = 1;
            } else if (s == "-null")  {
                shrinkposw = 0;
            } else if (s == "-shrinkposw")  {
                i++;
                shrinkposw = atof(argv[i]);
            } else if (s == "-shrinkdposom")  {
                i++;
                shrinkdposom = atof(argv[i]);
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

    CodonM2aSample *sample = new CodonM2aSample(newpath, name, burnin, every, until);
    if (ppred) {
        sample->PostPredSimu(shrinkposw,shrinkdposom);
    } else {
        sample->Read();
    }
}
