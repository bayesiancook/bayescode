#include <cmath>
#include <fstream>
#include <algorithm>
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
        vector<double> mean_branchsyn_array(Nbranch,0);
        vector<double> mean_branchom_array(Nbranch,0);
        vector<double> mean_genesyn_array(Ngene,0);
        vector<double> mean_geneom_array(Ngene,0);

        cerr << size << " points to read\n";
        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->AddBranchSynArrayTo(mean_branchsyn_array);
            GetModel()->AddBranchOmegaArrayTo(mean_branchom_array);
            GetModel()->AddGeneSynArrayTo(mean_genesyn_array);
            GetModel()->AddGeneOmegaArrayTo(mean_geneom_array);
        }
        cerr << '\n';
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
            GetModel()->AddSynDevLogFactorTo(syn_logfactor);
            GetModel()->AddOmegaDevLogFactorTo(om_logfactor);
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

        /*
        ofstream os((name + ".devpostprob").c_str());
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
        */

        ofstream os((name + ".devpostprob").c_str());
        os << "cutoff ratio and pp : " << exp(logratio_cutoff) << '\t' << pp_cutoff << '\n';
        for (int i=0; i<Ngene; i++)   {
            int nsyn = 0;
            int nom = 0;
            int nboth = 0;
            for (int j=0; j<Nbranch; j++)   {
                bool csyn = (fabs(syn_logfactor[i][j]) >= logratio_cutoff) && (syn_postprob[i][j] <= pp_cutoff);
                bool com = (fabs(om_logfactor[i][j]) >= logratio_cutoff) && (om_postprob[i][j] <= pp_cutoff);
                if (csyn)   {
                    nsyn++;
                }
                if (com)    {
                    nom++;
                }
                if (csyn && com)    {
                    nboth++;
                }
            }
            if (nsyn || nom)  {
                os << GetModel()->GetGeneName(i) << '\t' << nsyn << '\t' << nom << '\t' << nboth;
                for (int j=0; j<Nbranch; j++)   {
                    bool csyn = (fabs(syn_logfactor[i][j]) >= logratio_cutoff) && (syn_postprob[i][j] <= pp_cutoff);
                    bool com = (fabs(om_logfactor[i][j]) >= logratio_cutoff) && (om_postprob[i][j] <= pp_cutoff);
                    if (csyn || com)    {
                        os << '\t' << j;
                        if (csyn)   {
                            os << '\t' << double(int(10*exp(syn_logfactor[i][j])))/10;
                        }
                        else    {
                            os << '\t' << " - ";
                        }
                        if (com)    {
                            os << '\t' << double(int(10*exp(om_logfactor[i][j])))/10;
                        }
                        else    {
                            os << '\t' << " - ";
                        }
                    }
                }
                os << '\n';
            }
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
            GetModel()->AddSynDevToHist(syn_postsample, syn_predsample, i*Ngene*Nbranch);
            GetModel()->AddOmegaDevToHist(om_postsample, om_predsample, i*Ngene*Nbranch);
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

    FastGeneBranchOmegaSample *sample = new FastGeneBranchOmegaSample(name, burnin, every, until);
    if (qqplot) {
        sample->ReadQQPlot();
    }
    else if (ppdev) {
        sample->ReadMixDev(log(ratio), cutoff);
    }
    else    {
        sample->Read();
    }
}
