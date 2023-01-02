#include <cmath>
#include <fstream>
#include <algorithm>
#include "Sample.hpp"
#include "BranchNucModel.hpp"
// #include "RecursiveNewick.hpp"

using namespace std;

/**
 * \brief An MCMC sample for BranchNucModel
 *
 * implements a simple read function, returning the MCMC estimate of the
 * posterior mean and standard deviation of omega=dN/dS
 */

class GeneBranchStrandSymmetricSample : public Sample {
  private:
    string modeltype;
    string datafile;
    string treefile;
    string taxonfile;
    int syn_devmode;

  public:
    string GetModelType() override { return modeltype; }

    BranchNucModel *GetModel() override { return (BranchNucModel *)model; }

    //! \brief Constructor (file name, burn-in, thinning and upper limit, see
    //! Sample)
    GeneBranchStrandSymmetricSample(string filename, int inburnin, int inevery, int inuntil)
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
        is >> datafile >> treefile >> taxonfile;
        is >> syn_devmode;
        int check;
        is >> check;
        if (check)  {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "BRANCHNUC") {
            model = new BranchNucModel(datafile, treefile, taxonfile, syn_devmode);
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
        vector<double> mean_branchsyn_array(Nbranch,0);
        vector<vector<double>> mean_branchnuc_array(Nbranch,vector<double>(5,0));

        double meanAC, meanAG, meanCA, meanCG, meanCT;
        double varAC, varAG, varCA, varCG, varCT;
        meanAC = meanAG = meanCA = meanCG = meanCT = 0;
        varAC = varAG = varCA = varCG = varCT = 0;

        cerr << size << " points to read\n";
        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->GetSynModel()->AddBranchArrayTo(mean_branchsyn_array);
            GetModel()->AddNucStatsTo(meanAC, meanAG, meanCA, meanCG, meanCT,
                    varAC, varAG, varCA, varCG, varCT);
            GetModel()->AddNucBranchRatesTo(mean_branchnuc_array);
        }
        cerr << '\n';

        meanAC /= size;
        meanAG /= size;
        meanCA /= size;
        meanCG /= size;
        meanCT /= size;

        varAC /= size;
        varAG /= size;
        varCA /= size;
        varCG /= size;
        varCT /= size;

        ofstream nos((name + ".nucstats").c_str());
        nos << "          \tmean\trel-var\n";
        nos << fixed << setw(5) << setprecision(2);
        nos << "A:T -> C:G\t" << meanAC << '\t' << varAC << '\n';
        nos << "A:T -> G:C\t" << meanAG << '\t' << varAG << '\n';
        nos << "C:G -> A:T\t" << meanCA << '\t' << varCA << '\n';
        nos << "C:G -> G:C\t" << meanCG << '\t' << varCG << '\n';
        nos << "C:G -> T:A\t" << meanCT << '\t' << varCT << '\n';
        cerr << "post mean nuc stats in " << name << ".nucstats\n";

        for (int j=0; j<Nbranch; j++)   {
            for (int k=0; k<5; k++) {
                mean_branchnuc_array[j][k] /= size;
            }
        }

        auto nuc2string = [&nuc=mean_branchnuc_array] (int j)   {

            ostringstream s;
            s << fixed << setprecision(2);
            s << nuc[j][0] << "_" << nuc[j][1] << "_" << nuc[j][2] << "_" << nuc[j][3] << "_" << nuc[j][4];
            return s.str();
        };

        auto syn2string = [&syn=mean_branchsyn_array] (int j)   {
            ostringstream s;
            s << syn[j];
            return s.str();
        };

        ofstream nucos((name + ".postmeannuc.tre").c_str());
        ToNewick(nucos, *GetModel()->GetTree(), syn2string, nuc2string);
        cerr << "post mean tree in " << name << ".postmeannuc.tre\n";

        for (int j=0; j<Nbranch; j++)   {
            mean_branchsyn_array[j] /= size;
        }

        ofstream tos((name + ".postmeansyn.tre").c_str());
        ToNewick(tos, *GetModel()->GetTree(), mean_branchsyn_array, mean_branchsyn_array);
        cerr << "post mean tree in " << name << ".postmeansyn.tre\n";
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

    GeneBranchStrandSymmetricSample *sample = new GeneBranchStrandSymmetricSample(name, burnin, every, until);
    sample->Read();
}
