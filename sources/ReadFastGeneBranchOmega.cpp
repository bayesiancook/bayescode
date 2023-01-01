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

        vector<vector<double>> syn_z(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> om_z(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> syn_z2(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> om_z2(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> syn_zerr(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> om_zerr(Ngene, vector<double>(Nbranch,0));

        double mean_syn_invshape = 0;
        double mean_om_invshape = 0;

        cerr << size << " points to read\n";
        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();
            GetModel()->GetSynModel()->AddBranchArrayTo(mean_branchsyn_array);
            GetModel()->GetOmegaModel()->AddBranchArrayTo(mean_branchom_array);
            GetModel()->GetSynModel()->AddGeneArrayTo(mean_genesyn_array);
            GetModel()->GetOmegaModel()->AddGeneArrayTo(mean_geneom_array);
            mean_syn_invshape += GetModel()->GetSynModel()->GetDevInvShape();
            mean_om_invshape += GetModel()->GetOmegaModel()->GetDevInvShape();

            GetModel()->GetSynModel()->AddZscoreTo(syn_z);
            GetModel()->GetOmegaModel()->AddZscoreTo(om_z);
            GetModel()->GetSynModel()->AddRelVarTo(syn_z2);
            GetModel()->GetOmegaModel()->AddRelVarTo(om_z2);
        }
        cerr << '\n';

        mean_syn_invshape /= size;
        mean_om_invshape /= size;
        cerr << "syn rel var : " << mean_syn_invshape << '\n';
        cerr << "om rel var  : " << mean_om_invshape << '\n';
        cerr << '\n';

        for (int j=0; j<Nbranch; j++)   {
            mean_branchsyn_array[j] /= size;
            mean_branchom_array[j] /= size;
        }
        for (int i=0; i<Ngene; i++) {
            mean_genesyn_array[i] /= size;
            mean_geneom_array[i] /= size;
        }

        for (int i=0; i<Ngene; i++)   {
            for (int j=0; j<Nbranch; j++)   {
                syn_z[i][j] /= size;
                om_z[i][j] /= size;
                syn_z2[i][j] /= size;
                om_z2[i][j] /= size;
                syn_zerr[i][j] = sqrt(syn_z2[i][j] - syn_z[i][j]*syn_z[i][j]);
                om_zerr[i][j] = sqrt(om_z2[i][j] - om_z[i][j]*om_z[i][j]);
            }
        }

        vector<double> branch_syn_relvar(Nbranch, 0);
        vector<double> branch_om_relvar(Nbranch, 0);
        vector<double> gene_syn_relvar(Ngene, 0);
        vector<double> gene_om_relvar(Ngene, 0);
        for (int i=0; i<Ngene; i++)   {
            for (int j=0; j<Nbranch; j++)   {
                branch_syn_relvar[j] += syn_z2[i][j];
                gene_syn_relvar[i] += syn_z2[i][j];
                branch_om_relvar[j] += om_z2[i][j];
                gene_om_relvar[i] += om_z2[i][j];
            }
        }
        for (int i=0; i<Ngene; i++) {
            gene_syn_relvar[i] /= Nbranch;
            gene_om_relvar[i] /= Nbranch;
        }
        for (int j=0; j<Nbranch; j++)   {
            branch_syn_relvar[j] /= Ngene;
            branch_om_relvar[j] /= Ngene;
        }

        ofstream devos((name + ".postmeandev.tab").c_str());
        devos << "#branchsynmean\tbranchsynrelvar\tgenesynmean\tgenesynrelvar\tsynz\tsynzerr\tbranchommean\tbranchomrelvar\tgeneommean\tgeneomrelvar\tomz\tomzerr\n";
        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                devos << mean_branchsyn_array[j] << '\t' << branch_syn_relvar[j] << '\t' << mean_genesyn_array[i] << '\t' << gene_syn_relvar[i] << '\t' << syn_z[i][j] << '\t' << syn_zerr[i][j] << '\t' << mean_branchom_array[j] << '\t' << branch_om_relvar[j] << '\t' << mean_geneom_array[i] << '\t' << gene_om_relvar[i] << '\t' << om_z[i][j] << '\t' << om_zerr[i][j] << '\n';
            }
        }

        ofstream gos((name + ".postmean.genesynom.tab").c_str());
        for (int i=0; i<Ngene; i++) {
            gos << GetModel()->GetGeneName(i) << '\t' << mean_genesyn_array[i] << '\t' << gene_syn_relvar[i] << '\t' << mean_geneom_array[i] << '\t' << gene_om_relvar[i] << '\n';
        }
        cerr << "post mean gene dN/dS in " << name << ".postmean.geneom.tab\n";

        ofstream bsos((name + ".postmean.branchsyn.tab").c_str());
        Tabulate(bsos, GetModel()->GetTree(), mean_branchsyn_array, branch_syn_relvar, false);

        ofstream boos((name + ".postmean.branchom.tab").c_str());
        Tabulate(boos, GetModel()->GetTree(), mean_branchom_array, branch_om_relvar, false);

        ofstream os((name + ".postmean.leafdsom.tab").c_str());
        Tabulate(os, GetModel()->GetTree(), mean_branchsyn_array, mean_branchom_array, true);

        ofstream tos((name + ".dsom.tre").c_str());
        ToNewick(tos, *GetModel()->GetTree(), mean_branchsyn_array, mean_branchom_array);

        cerr << "post mean branch effects on dN/dS in " << name << ".postmean.leafdsom.tab\n";
        cerr << "newick format in " << name << ".dsom.tre\n";
        cerr << '\n';
    }

    void ReadMixDev(double z_cutoff, double pp_cutoff)   {
        int Ngene = GetModel()->GetNgene();
        int Nbranch = GetModel()->GetNbranch();
        vector<vector<double>> syn_postprob(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> om_postprob(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> syn_z(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> om_z(Ngene, vector<double>(Nbranch,0));

        double syn_devpi = 0;
        double om_devpi = 0;

        cerr << size << " points to read\n";
        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();
            GetModel()->AddSynDevPostProbsTo(syn_postprob);
            GetModel()->AddOmegaDevPostProbsTo(om_postprob);
            GetModel()->GetSynModel()->AddDevZscoreTo(syn_z);
            GetModel()->GetOmegaModel()->AddDevZscoreTo(om_z);

            syn_devpi += GetModel()->GetSynModel()->GetDevPi();
            om_devpi += GetModel()->GetOmegaModel()->GetDevPi();
        }
        cerr << '\n';

        syn_devpi /= size;
        om_devpi /= size;

        for (int i=0; i<Ngene; i++)   {
            for (int j=0; j<Nbranch; j++)   {
                syn_postprob[i][j] /= size;
                om_postprob[i][j] /= size;
                syn_z[i][j] /= size;
                om_z[i][j] /= size;
            }
        }

        ofstream os((name + ".devzscores").c_str());
        int totsyn = 0;
        int totom = 0;
        int totgenesyn = 0;
        int totgeneom = 0;
        for (int i=0; i<Ngene; i++)   {
            os << GetModel()->GetGeneName(i);
            int syn_one = 0;
            int om_one = 0;
            for (int j=0; j<Nbranch; j++)   {
                if ((fabs(syn_z[i][j] > z_cutoff)) && (syn_postprob[i][j] < pp_cutoff))    {
                    totsyn++;
                    syn_one = 1;
                }
                if ((fabs(om_z[i][j] > z_cutoff)) && (om_postprob[i][j] < pp_cutoff))    {
                    totom++;
                    om_one = 1;
                }
                os << '\t' << syn_postprob[i][j];
                os << '\t' << syn_z[i][j];
                os << '\t' << om_postprob[i][j];
                os << '\t' << om_z[i][j];
            }
            os << '\n';
            if (syn_one)    {
                totgenesyn++;
            }
            if (om_one) {
                totgeneom++;
            }
        }
        int exptotsyn = syn_devpi * Nbranch * Ngene;
        int exptotom = om_devpi * Nbranch * Ngene;
        cerr << "number of deviating gene/branch effects\n";
        cerr << "syn : " << totsyn << " (" << exptotsyn << ")" << '\t' << totgenesyn << '\n';
        cerr << "om  : " << totom << " (" << exptotom << ") " << '\t' << totgeneom << '\n';
    }

    void ReadQQPlot(double z_cutoff = 1.0)   {
        int Ngene = GetModel()->GetNgene();
        int Nbranch = GetModel()->GetNbranch();
        int n = Ngene*Nbranch*size;
        vector<double> syn_postsample(n,0);
        vector<double> syn_predsample(n,0);
        vector<double> om_postsample(n,0);
        vector<double> om_predsample(n,0);

        vector<vector<double>> syn_z(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> om_z(Ngene, vector<double>(Nbranch,0));

        cerr << size << " points to read\n";
        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();
            GetModel()->GetSynModel()->AddDevToHist(syn_postsample, syn_predsample, i*Ngene*Nbranch);
            GetModel()->GetOmegaModel()->AddDevToHist(om_postsample, om_predsample, i*Ngene*Nbranch);
            GetModel()->GetSynModel()->AddDevZscoreTo(syn_z);
            GetModel()->GetOmegaModel()->AddDevZscoreTo(om_z);
        }
        cerr << '\n';
        ofstream os((name + ".qqplot").c_str());
        sort(syn_postsample.begin(), syn_postsample.end());
        sort(syn_predsample.begin(), syn_predsample.end());
        sort(om_postsample.begin(), om_postsample.end());
        sort(om_predsample.begin(), om_predsample.end());
        int synposti = 0;
        int synpredi = 0;
        int omposti = 0;
        int ompredi = 0;
        for (int i=0; i<n; i++) {
            os << syn_postsample[i] << '\t' << syn_predsample[i] << '\t' << om_postsample[i] << '\t' << om_predsample[i] << '\n';
            if (syn_postsample[i] < z_cutoff)   {
                synposti = i;
            }
            if (syn_predsample[i] < z_cutoff)   {
                synpredi = i;
            }
            if (om_postsample[i] < z_cutoff)    {
                omposti = i;
            }
            if (om_predsample[i] < z_cutoff)    {
                ompredi = i;
            }
        }
        cerr << "sorted centered and normalized deviations in " << name << ".qqplot\n";
        int synpostf = (n-synposti-1);
        int synpredf = (n-synpredi-1);
        int ompostf = (n-omposti-1);
        int ompredf = (n-ompredi-1);
        cerr << "fraction with z > " << z_cutoff << " (obs / exp)\n";
        cerr << "syn : " << synpostf << " / " << synpredf << '\n';
        cerr << "om  : " << ompostf << " / " << ompredf << '\n';

        for (int i=0; i<Ngene; i++)   {
            for (int j=0; j<Nbranch; j++)   {
                syn_z[i][j] /= size;
                om_z[i][j] /= size;
            }
        }
        ofstream los((name + ".devzscores").c_str());
        int totsyn = 0;
        int totom = 0;
        int totgenesyn = 0;
        int totgeneom = 0;
        for (int i=0; i<Ngene; i++)   {
            los << GetModel()->GetGeneName(i);
            int syn_one = 0;
            int om_one = 0;
            for (int j=0; j<Nbranch; j++)   {
                if (syn_z[i][j] > z_cutoff)  {
                    totsyn++;
                    syn_one = 1;
                }
                if (om_z[i][j] > z_cutoff)   {
                    totom++;
                    om_one = 1;
                }
                los << '\t' << syn_z[i][j];
                los << '\t' << om_z[i][j];
            }
            los << '\n';
            if (syn_one)    {
                totgenesyn++;
            }
            if (om_one) {
                totgeneom++;
            }
        }
        cerr << "number of deviating gene/branch effects\n";
        cerr << "syn : " << totsyn << '\t' << totgenesyn << '\n';
        cerr << "om  : " << totom << '\t' << totgeneom << '\n';
    }
};

int main(int argc, char *argv[]) {
    int burnin = 0;
    int every = 1;
    int until = -1;
    int dev = 0;
    int ppdev = 0;
    double z = 0;
    double cutoff = 0;

    string name;

    try {
        if (argc == 1) {
            throw(0);
        }

        int i = 1;
        while (i < argc) {
            string s = argv[i];
            if (s == "-dev") {
                dev = 1;
                i++;
                z = atof(argv[i]);
            } else if (s == "-mixdev")   {
                ppdev = 1;
                i++;
                z = atof(argv[i]);
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
    if (dev) {
        sample->ReadQQPlot(z);
    }
    else if (ppdev) {
        sample->ReadMixDev(z, cutoff);
    }
    else    {
        sample->Read();
    }
}
