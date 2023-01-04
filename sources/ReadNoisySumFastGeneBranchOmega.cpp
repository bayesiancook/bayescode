#include <cmath>
#include <fstream>
#include <algorithm>
#include "Sample.hpp"
#include "NoisySumFastGeneBranchOmegaModel.hpp"
// #include "RecursiveNewick.hpp"

using namespace std;

/**
 * \brief An MCMC sample for NoisySumFastGeneBranchOmegaModel
 *
 * implements a simple read function, returning the MCMC estimate of the
 * posterior mean and standard deviation of omega=dN/dS
 */

class NoisySumFastGeneBranchOmegaSample : public Sample {
  private:
    string modeltype;
    string datafile;
    string treefile;

  public:
    string GetModelType() override { return modeltype; }

    NoisySumFastGeneBranchOmegaModel *GetModel() override { return (NoisySumFastGeneBranchOmegaModel *)model; }

    //! \brief Constructor (file name, burn-in, thinning and upper limit, see
    //! Sample)
    NoisySumFastGeneBranchOmegaSample(string filename, int inburnin, int inevery, int inuntil)
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
        int check;
        is >> check;
        if (check)  {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "NOISYSUMFASTGENEBRANCHOMEGA") {
            model = new NoisySumFastGeneBranchOmegaModel(datafile, treefile);
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
    void Read(double z_cutoff, double pp_cutoff)   {

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

        vector<vector<double>> syn_postprob(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> om_postprob(Ngene, vector<double>(Nbranch,0));

        double syn_invshape = 0;
        double syn_mean2 = 0;
        double syn_invshape2 = 0;

        double om_invshape = 0;
        double om_mean2 = 0;
        double om_invshape2 = 0;

        cerr << size << " points to read\n";
        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();

            GetModel()->GetSynModel()->AddBranchArrayTo(mean_branchsyn_array);
            GetModel()->GetOmegaModel()->AddBranchArrayTo(mean_branchom_array);
            GetModel()->GetSynModel()->AddGeneArrayTo(mean_genesyn_array);
            GetModel()->GetOmegaModel()->AddGeneArrayTo(mean_geneom_array);

            syn_invshape += GetModel()->GetSynModel()->GetDevInvShape();
            syn_mean2 += GetModel()->GetSynModel()->GetDevMean2();
            syn_invshape2 += GetModel()->GetSynModel()->GetDevInvShape2();

            om_invshape += GetModel()->GetOmegaModel()->GetDevInvShape();
            om_mean2 += GetModel()->GetOmegaModel()->GetDevMean2();
            om_invshape2 += GetModel()->GetOmegaModel()->GetDevInvShape2();

            GetModel()->GetSynModel()->AddZscoreTo(syn_z);
            GetModel()->GetOmegaModel()->AddZscoreTo(om_z);

            GetModel()->GetSynModel()->AddRelVarTo(syn_z2);
            GetModel()->GetOmegaModel()->AddRelVarTo(om_z2);

            GetModel()->AddSynDevPostProbsTo(syn_postprob);
            GetModel()->AddOmegaDevPostProbsTo(om_postprob);
        }
        cerr << '\n';

        syn_invshape /= size;
        syn_mean2 /= size;
        syn_invshape2 /= size;

        om_invshape /= size;
        om_mean2 /= size;
        om_invshape2 /= size;

        cerr << "\trelvar\tmean2\trelvar2\n";
        cerr << "syn" << syn_invshape << '\t' << syn_mean2 << '\t' << syn_invshape2 << '\n';
        cerr << "om" << om_invshape << '\t' << om_mean2 << '\t' << om_invshape2 << '\n';
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
                syn_postprob[i][j] /= size;
                om_postprob[i][j] /= size;

                syn_z[i][j] /= size;
                om_z[i][j] /= size;

                syn_z2[i][j] /= size;
                om_z2[i][j] /= size;
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
        devos << "#genename\tbranchsynmean\tbranchsynrelvar\tgenesynmean\tgenesynrelvar\tsynz\tsynpp\tbranchommean\tbranchomrelvar\tgeneommean\tgeneomrelvar\tomz\tompp\n";

        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                devos << GetModel()->GetGeneName(i);
                devos << '\t' << mean_branchsyn_array[j] << '\t' << branch_syn_relvar[j];
                devos << '\t' << mean_genesyn_array[i] << '\t' << gene_syn_relvar[i];
                devos << '\t' << syn_z[i][j] << '\t' << syn_postprob[i][j];
                devos << '\t' << mean_branchom_array[j] << '\t' << branch_om_relvar[j];
                devos << '\t' << mean_geneom_array[i] << '\t' << gene_om_relvar[i];
                devos << '\t' << om_z[i][j] << '\t' << om_postprob[i][j];
                devos << '\n';
            }
        }

        ofstream gos((name + ".postmean.genesynom.tab").c_str());
        gos << "#genename\tgenesynmean\tgenesynrelvar\tgeneommean\tgeneomrelvar\n";
        for (int i=0; i<Ngene; i++) {
            gos << GetModel()->GetGeneName(i) << '\t' << mean_genesyn_array[i] << '\t' << gene_syn_relvar[i] << '\t' << mean_geneom_array[i] << '\t' << gene_om_relvar[i] << '\n';
        }
        cerr << "post mean gene dN/dS in " << name << ".postmean.geneom.tab\n";

        ofstream bsos((name + ".postmean.branchsyn.tab").c_str());
        bsos << "#taxon1\ttaxon2\tbranchsynmean\tbranchsynrelvar\n";
        Tabulate(bsos, GetModel()->GetTree(), mean_branchsyn_array, branch_syn_relvar, false);

        ofstream boos((name + ".postmean.branchom.tab").c_str());
        boos << "#taxon1\ttaxon2\tbranchommean\tbranchomrelvar\n";
        Tabulate(boos, GetModel()->GetTree(), mean_branchom_array, branch_om_relvar, false);

        ofstream os((name + ".postmean.leafdsom.tab").c_str());
        os << "#taxon\tbranchsyn\tbranchom\n";
        Tabulate(os, GetModel()->GetTree(), mean_branchsyn_array, mean_branchom_array, true);

        ofstream tos((name + ".dsom.tre").c_str());
        ToNewick(tos, *GetModel()->GetTree(), mean_branchsyn_array, mean_branchom_array);

        cerr << "post mean branch effects on dN/dS in " << name << ".postmean.leafdsom.tab\n";
        cerr << "newick format in " << name << ".dsom.tre\n";
        cerr << '\n';

        int totsyn = 0;
        int totom = 0;
        vector<int> syn_genecount(Ngene,0);
        vector<int> om_genecount(Ngene,0);
        for (int i=0; i<Ngene; i++)   {
            for (int j=0; j<Nbranch; j++)   {
                if ((fabs(syn_z[i][j] > z_cutoff)) && (syn_postprob[i][j] < pp_cutoff))    {
                    totsyn++;
                    syn_genecount[i]++;
                }
                if ((fabs(om_z[i][j] > z_cutoff)) && (om_postprob[i][j] < pp_cutoff))    {
                    totom++;
                    om_genecount[i]++;
                }
            }
        }
        /*
        int exptotsyn = syn_pi * Nbranch * Ngene;
        int exptotom = om_pi * Nbranch * Ngene;
        int totgenesyn = 0;
        int totgeneom = 0;
        int totgeneboth = 0;
        for (int i=0; i<Ngene; i++)	{
            if (syn_genecount[i])	{
                totgenesyn++;
            }
            if (om_genecount[i])    {
                totgeneom++;
            }
            if (om_genecount[i] || syn_genecount[i])    {
                totgeneboth++;
            }
        }

        cerr << "number of deviating gene/branch effects\n";
        cerr << "syn : " << totsyn << " (" << exptotsyn << ")" << '\t' << totgenesyn << '\n';
        cerr << "om  : " << totom << " (" << exptotom << ") " << '\t' << totgeneom << '\n';
        cerr << "both: " << totgeneboth << '\n';
        */
    }
};

int main(int argc, char *argv[]) {
    int burnin = 0;
    int every = 1;
    int until = -1;
    double z = 5;
    double cutoff = 0.05;

    string name;

    try {
        if (argc == 1) {
            throw(0);
        }

        int i = 1;
        while (i < argc) {
            string s = argv[i];
            if (s == "-z")   {
                i++;
                z = atof(argv[i]);
            } else if (s == "-pp")	{
                i++;
                cutoff = atof(argv[i]);
            } else if ((s == "-x") || (s == "-extract")) {
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

    NoisySumFastGeneBranchOmegaSample *sample = new NoisySumFastGeneBranchOmegaSample(name, burnin, every, until);
    sample->Read(z, cutoff);
}
