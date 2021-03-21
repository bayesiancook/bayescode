
#include <cmath>
#include <fstream>
#include "Sample.hpp"
#include "FastCoevolModel.hpp"
#include "DistBranchNodeArray.hpp"
#include "MeanCovMatrix.hpp"
using namespace std;

/**
 * \brief An MCMC sample for FastCoevolModel
 *
 * implements a simple read function, returning the MCMC estimate of the
 * posterior mean and standard deviation of omega=dN/dS
 */

class FastCoevolSample : public Sample {
  private:
    string modeltype;
    string contdatafile, treefile, rootfile;
    string dsomsuffstatfile;
    int wndsmode, wnommode;

  public:
    string GetModelType() override { return modeltype; }

    FastCoevolModel *GetModel() override { return (FastCoevolModel *)model; }

    //! \brief Constructor (file name, burn-in, thinning and upper limit, see
    //! Sample)
    FastCoevolSample(string filename, int inburnin, int inevery, int inuntil)
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
        is >> contdatafile >> treefile >> rootfile;
        is >> dsomsuffstatfile;
        is >> wndsmode >> wnommode;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "FASTCOEVOLDNDS") {
            model = new FastCoevolModel(contdatafile, treefile, rootfile, dsomsuffstatfile, wndsmode, wnommode);
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

        DistBranchNodeArray meansynrate(GetModel()->GetTree());
        DistBranchNodeArray meanomega(GetModel()->GetTree());

		int dim = GetModel()->GetCovMatrix().GetDim();
		MeanCovMatrix  mat(dim);

        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();
            meansynrate.AddFromChrono(GetModel()->GetChronogram(), GetModel()->GetProcess(), 0);
            meanomega.AddFromChrono(GetModel()->GetChronogram(), GetModel()->GetProcess(), 1);
			mat.Add(GetModel()->GetCovMatrix());
        }
        cerr << '\n';
    
        meansynrate.Sort();
        ofstream sos((name + ".postmeands.tre").c_str());
        meansynrate.MedianToStream(sos);
        cerr << "postmean dS tree in " << name << ".postmeands.tre\n"; 

        meanomega.Sort();
        ofstream omos((name + ".postmeanomega.tre").c_str());
        meanomega.MedianToStream(omos);
        cerr << "postmean omega tree in " << name << ".postmeanomega.tre\n"; 

		mat.Normalize();
		ofstream mos((name + ".cov").c_str());
		mos << "entries are in the following order:\n";
		GetModel()->PrintEntries(mos);
		mos << '\n';
		// mat.SetLatex(tex);
		mos << mat;
		cerr << "covariance matrix in " << name << ".cov\n";
		cerr << '\n';
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

    FastCoevolSample *sample = new FastCoevolSample(name, burnin, every, until);
    if (ppred) {
        sample->PostPred();
    } else {
        sample->Read();
    }
}
