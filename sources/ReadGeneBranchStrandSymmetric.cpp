#include <cmath>
#include <fstream>
#include <algorithm>
#include "Sample.hpp"
#include "GeneBranchStrandSymmetricCodonModel.hpp"
// #include "RecursiveNewick.hpp"

using namespace std;

/**
 * \brief An MCMC sample for GeneBranchStrandSymmetricCodonModel
 *
 * implements a simple read function, returning the MCMC estimate of the
 * posterior mean and standard deviation of omega=dN/dS
 */

class GeneBranchStrandSymmetricSample : public Sample {
  private:
    string modeltype;
    string datafile;
    string treefile;
    int syn_devmode, om_devmode, nuc_devmode;

  public:
    string GetModelType() override { return modeltype; }

    GeneBranchStrandSymmetricCodonModel *GetModel() override { return (GeneBranchStrandSymmetricCodonModel *)model; }

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
        is >> datafile >> treefile;
        is >> syn_devmode >> om_devmode >> nuc_devmode;
        int check;
        is >> check;
        if (check)  {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "GENEBRANCHCODON") {
            model = new GeneBranchStrandSymmetricCodonModel(datafile, treefile, syn_devmode, om_devmode, nuc_devmode);
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
        vector<double> mean_genesyn_array(Ngene,0);
        vector<double> mean_geneom_array(Ngene,0);

        double meanAC, meanAG, meanCA, meanCG, meanCT;
        double geneAC, geneAG, geneCA, geneCG, geneCT;
        double branchAC, branchAG, branchCA, branchCG, branchCT;
        double devAC, devAG, devCA, devCG, devCT;

        meanAC = meanAG = meanCA = meanCG = meanCT = 0;
        geneAC = geneAG = geneCA = geneCG = geneCT = 0;
        branchAC = branchAG = branchCA = branchCG = branchCT = 0;
        devAC = devAG = devCA = devCG = devCT = 0;

        cerr << size << " points to read\n";
        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->GetSynModel()->AddBranchArrayTo(mean_branchsyn_array);
            GetModel()->GetOmegaModel()->AddBranchArrayTo(mean_branchom_array);
            GetModel()->GetSynModel()->AddGeneArrayTo(mean_genesyn_array);
            GetModel()->GetOmegaModel()->AddGeneArrayTo(mean_geneom_array);
            GetModel()->AddNucStats(meanAC, meanAG, meanCA, meanCG, meanCT,
                    geneAC, geneAG, geneCA, geneCG, geneCT,
                    branchAC, branchAG, branchCA, branchCG, branchCT,
                    devAC, devAG, devCA, devCG, devCT);
        }
        cerr << '\n';

        meanAC /= size;
        meanAG /= size;
        meanCA /= size;
        meanCG /= size;
        meanCT /= size;

        geneAC /= size;
        geneAG /= size;
        geneCA /= size;
        geneCG /= size;
        geneCT /= size;

        branchAC /= size;
        branchAG /= size;
        branchCA /= size;
        branchCG /= size;
        branchCT /= size;

        devAC /= size;
        devAG /= size;
        devCA /= size;
        devCG /= size;
        devCT /= size;

        ofstream nos((name + "nucstats").c_str());
        nos << "XX\tmean\tgenevar\tbranchvar\tdevvar\n";
        nos << meanAC << '\t' << geneAC << '\t' << branchAC << '\t' << devAC << '\n';
        nos << meanAG << '\t' << geneAG << '\t' << branchAG << '\t' << devAG << '\n';
        nos << meanCA << '\t' << geneCA << '\t' << branchCA << '\t' << devCA << '\n';
        nos << meanCG << '\t' << geneCG << '\t' << branchCG << '\t' << devCG << '\n';
        nos << meanCT << '\t' << geneCT << '\t' << branchCT << '\t' << devCT << '\n';
        cerr << "post mean nuc stats in " << name << ".nuctstats\n";

        for (int j=0; j<Nbranch; j++)   {
            mean_branchsyn_array[j] /= size;
            mean_branchom_array[j] /= size;
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

        ofstream tos((name + ".dsom.tre").c_str());
        ToNewick(tos, *GetModel()->GetTree(), mean_branchsyn_array, mean_branchom_array);

        cerr << "post mean branch effects on dN/dS in " << name << ".postmean.leafdsom.tab\n";
        cerr << "newick format in " << name << ".dsom.tre\n";
    }

    /*
    void ReadMixDev(double logratio_cutoff, double pp_cutoff)   {
        int Ngene = GetModel()->GetNgene();
        int Nbranch = GetModel()->GetNbranch();
        vector<vector<double>> syn_postprob(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> om_postprob(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> syn_logfactor(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> om_logfactor(Ngene, vector<double>(Nbranch,0));
        cerr << size << " points to read\n";
        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();
            GetModel()->AddSynDevPostProbsTo(syn_postprob);
            GetModel()->AddOmegaDevPostProbsTo(om_postprob);
            GetModel()->GetSynModel()->AddDevLogFactorTo(syn_logfactor);
            GetModel()->GetOmegaModel()->AddDevLogFactorTo(om_logfactor);
        }
        cerr << '\n';
        for (int i=0; i<Ngene; i++)   {
            for (int j=0; j<Nbranch; j++)   {
                syn_postprob[i][j] /= size;
                om_postprob[i][j] /= size;
                syn_logfactor[i][j] /= size;
                om_logfactor[i][j] /= size;
            }
        }

        ofstream os((name + ".devlogratios").c_str());
        for (int i=0; i<Ngene; i++)   {
            os << GetModel()->GetGeneName(i);
            for (int j=0; j<Nbranch; j++)   {
                os << '\t' << syn_postprob[i][j];
                os << '\t' << syn_logfactor[i][j];
                os << '\t' << om_postprob[i][j];
                os << '\t' << om_logfactor[i][j];
            }
            os << '\n';
        }
    }

    void ReadQQPlot()   {
        int Ngene = GetModel()->GetNgene();
        int Nbranch = GetModel()->GetNbranch();
        int n = Ngene*Nbranch*size;
        vector<double> syn_postsample(n,0);
        vector<double> syn_predsample(n,0);
        vector<double> om_postsample(n,0);
        vector<double> om_predsample(n,0);
        cerr << size << " points to read\n";
        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();
            GetModel()->GetSynModel()->AddDevToHist(syn_postsample, syn_predsample, i*Ngene*Nbranch);
            GetModel()->GetOmegaModel()->AddDevToHist(om_postsample, om_predsample, i*Ngene*Nbranch);
        }
        cerr << '\n';
        ofstream os((name + ".qqplot").c_str());
        sort(syn_postsample.begin(), syn_postsample.end());
        sort(syn_predsample.begin(), syn_predsample.end());
        sort(om_postsample.begin(), om_postsample.end());
        sort(om_predsample.begin(), om_predsample.end());
        for (int i=0; i<n; i++) {
            os << syn_postsample[i] << '\t' << syn_predsample[i] << '\t' << om_postsample[i] << '\t' << om_predsample[i] << '\n';
        }
        cerr << "sorted centered deviations in " << name << ".qqplot\n";
    }
    */
};

int main(int argc, char *argv[]) {
    int burnin = 0;
    int every = 1;
    int until = -1;
    int qqplot = 0;
    int ppdev = 0;
    double ratio = 0;
    double cutoff = 0;

    string name;

    try {
        if (argc == 1) {
            throw(0);
        }

        int i = 1;
        while (i < argc) {
            string s = argv[i];
            if (s == "-qqplot") {
                qqplot = 1;
            } else if (s == "-ppdev")   {
                ppdev = 1;
                i++;
                ratio = atof(argv[i]);
                i++;
                cutoff = atof(argv[i]);
            }
            else if ((s == "-x") || (s == "-extract")) {
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
    /*
    if (qqplot) {
        sample->ReadQQPlot();
    }
    else if (ppdev) {
        sample->ReadMixDev(log(ratio), cutoff);
    }
    else    {
        sample->Read();
    }
    */
}
