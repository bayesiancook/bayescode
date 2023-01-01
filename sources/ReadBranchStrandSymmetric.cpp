#include <cmath>
#include <fstream>
#include <algorithm>
#include "Sample.hpp"
#include "BranchStrandSymmetricCodonModel.hpp"
// #include "RecursiveNewick.hpp"

using namespace std;

/**
 * \brief An MCMC sample for BranchStrandSymmetricCodonModel
 *
 * implements a simple read function, returning the MCMC estimate of the
 * posterior mean and standard deviation of omega=dN/dS
 */

class BranchStrandSymmetricSample : public Sample {
  private:
    string modeltype;
    string datafile;
    string treefile;
    string taxonfile;
    int syn_devmode, om_devmode;

  public:
    string GetModelType() override { return modeltype; }

    BranchStrandSymmetricCodonModel *GetModel() override { return (BranchStrandSymmetricCodonModel *)model; }

    //! \brief Constructor (file name, burn-in, thinning and upper limit, see
    //! Sample)
    BranchStrandSymmetricSample(string filename, int inburnin, int inevery, int inuntil)
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
        is >> syn_devmode >> om_devmode;
        int check;
        is >> check;
        if (check)  {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "BRANCHCODON") {
            model = new BranchStrandSymmetricCodonModel(datafile, treefile, taxonfile, syn_devmode, om_devmode);
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
        vector<double> mean_branchsyn_array(Nbranch,0);
        vector<double> mean_branchom_array(Nbranch,0);
        vector<double> var_branchom_array(Nbranch,0);
        vector<double> err_branchom_array(Nbranch,0);
        vector<double> mean_genesyn_array(Ngene,0);
        vector<double> mean_geneom_array(Ngene,0);

        double meanAC, meanAG, meanCA, meanCG, meanCT;
        double varAC, varAG, varCA, varCG, varCT;
        meanAC = meanAG = meanCA = meanCG = meanCT = 0;
        varAC = varAG = varCA = varCG = varCT = 0;

        cerr << size << " points to read\n";
        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->GetSynModel()->AddBranchArrayTo(mean_branchsyn_array);
            GetModel()->GetOmegaModel()->AddBranchArrayTo(mean_branchom_array);
            GetModel()->GetOmegaModel()->AddSquaredBranchArrayTo(var_branchom_array);
            GetModel()->GetSynModel()->AddGeneArrayTo(mean_genesyn_array);
            GetModel()->GetOmegaModel()->AddGeneArrayTo(mean_geneom_array);
            GetModel()->AddNucStatsTo(meanAC, meanAG, meanCA, meanCG, meanCT,
                    varAC, varAG, varCA, varCG, varCT);
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
            mean_branchsyn_array[j] /= size;
            mean_branchom_array[j] /= size;
            var_branchom_array[j] /= size;
            var_branchom_array[j] -= mean_branchom_array[j] * mean_branchom_array[j];
            err_branchom_array[j] = sqrt(var_branchom_array[j]);
        }
        for (int i=0; i<Ngene; i++) {
            mean_genesyn_array[i] /= size;
            mean_geneom_array[i] /= size;
        }
        ofstream gos((name + ".postmean.geneom.tab").c_str());
        for (int i=0; i<Ngene; i++) {
            gos << GetModel()->GetGeneName(i) << '\t' << mean_genesyn_array[i] << '\t' << mean_geneom_array[i] << '\n';
        }
        cerr << "post mean gene dN/dS in " << name << ".postmean.geneom.tab\n";

        ofstream os((name + ".postmean.leafdsom.tab").c_str());
        Tabulate(os, GetModel()->GetTree(), mean_branchsyn_array, mean_branchom_array, true);
        // Tabulate(os, GetModel()->GetTree(), mean_branchom_array, err_branchom_array, true);

        ofstream tos((name + ".dsom.tre").c_str());
        ToNewick(tos, *GetModel()->GetTree(), mean_branchsyn_array, mean_branchom_array);

        cerr << "post mean branch effects on dN/dS in " << name << ".postmean.leafdsom.tab\n";
        cerr << "newick format in " << name << ".dsom.tre\n";
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

    BranchStrandSymmetricSample *sample = new BranchStrandSymmetricSample(name, burnin, every, until);
    sample->Read();
}
