#include <cmath>
#include <fstream>
#include "Sample.hpp"
#include "FastGeneBranchOmegaModel.hpp"
// #include "RecursiveNewick.hpp"

using namespace std;

/**
 * \brief An MCMC sample for FastGeneBranchOmegaModel
 *
 * implements a simple read function, returning the MCMC estimate of the
 * posterior mean and standard deviation of omega=dN/dS
 */

class FastGeneBranchOmegaSample : public Sample {
  private:
    string modeltype;
    string datafile;
    string treefile;
    int syn_devmode, om_devmode;

  public:
    string GetModelType() override { return modeltype; }

    FastGeneBranchOmegaModel *GetModel() override { return (FastGeneBranchOmegaModel *)model; }

    //! \brief Constructor (file name, burn-in, thinning and upper limit, see
    //! Sample)
    FastGeneBranchOmegaSample(string filename, int inburnin, int inevery, int inuntil)
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
        is >> syn_devmode >> om_devmode;
        int check;
        is >> check;
        if (check)  {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "FASTGENEBRANCHOMEGA") {
            model = new FastGeneBranchOmegaModel(datafile, treefile, syn_devmode, om_devmode);
        } else {
            cerr << "error when opening file " << name << '\n';
            exit(1);
        }

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
        int Nbranch = GetModel()->GetNbranch();
        int Ngene = GetModel()->GetNgene();
        vector<double> mean_branchom_array(Nbranch,0);
        vector<double> mean_geneom_array(Ngene,0);
        cerr << size << " points to read\n";
        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->AddBranchOmegaArrayTo(mean_branchom_array);
            GetModel()->AddGeneOmegaArrayTo(mean_geneom_array);
        }
        cerr << '\n';
        for (int j=0; j<Nbranch; j++)   {
            mean_branchom_array[j] /= size;
        }
        for (int i=0; i<Ngene; i++) {
            mean_geneom_array[i] /= size;
        }
        ofstream os((name + ".postmean.leafdsom.tab").c_str());
        Tabulate(os, GetModel()->GetTree(), mean_branchom_array, true);
        ofstream gos((name + ".postmean.geneom.tab").c_str());
        for (int i=0; i<Ngene; i++) {
            gos << GetModel()->GetGeneName(i) << '\t' << mean_geneom_array[i] << '\n';
        }
        cerr << "post mean gene dN/dS in " << name << ".postmean.geneom.tab\n";
        cerr << "post mean branch effects on dN/dS in " << name << ".postmean.leafdsom.tab\n";
    }
};

int main(int argc, char *argv[]) {
    int burnin = 0;
    int every = 1;
    int until = -1;

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

    FastGeneBranchOmegaSample *sample = new FastGeneBranchOmegaSample(name, burnin, every, until);
    sample->Read();
}
