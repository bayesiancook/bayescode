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

        vector<vector<double>> syn_val(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> om_val(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> syn_z(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> om_z(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> syn_z2(Ngene, vector<double>(Nbranch,0));
        vector<vector<double>> om_z2(Ngene, vector<double>(Nbranch,0));

        double mean_syn = 0;
        double gene_syn = 0;
        double branch_syn = 0;
        double dev_syn = 0;
        double mean_om = 0;
        double gene_om = 0;
        double branch_om = 0;
        double dev_om = 0;

        cerr << size << " points to read\n";
        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();
            GetModel()->GetSynModel()->AddBranchArrayTo(mean_branchsyn_array);
            GetModel()->GetOmegaModel()->AddBranchArrayTo(mean_branchom_array);
            GetModel()->GetSynModel()->AddGeneArrayTo(mean_genesyn_array);
            GetModel()->GetOmegaModel()->AddGeneArrayTo(mean_geneom_array);

            GetModel()->GetSynModel()->AddValTo(syn_val);
            GetModel()->GetOmegaModel()->AddValTo(om_val);
            GetModel()->GetSynModel()->AddZscoreTo(syn_z);
            GetModel()->GetOmegaModel()->AddZscoreTo(om_z);
            GetModel()->GetSynModel()->AddRelVarTo(syn_z2);
            GetModel()->GetOmegaModel()->AddRelVarTo(om_z2);

            GetModel()->GetSynModel()->AddStats(mean_syn, gene_syn, branch_syn, dev_syn);
            GetModel()->GetOmegaModel()->AddStats(mean_om, gene_om, branch_om, dev_om);
        }
        cerr << '\n';

        mean_syn /= size;
        gene_syn /= size;
        branch_syn /= size;
        dev_syn /= size;
        mean_om /= size;
        gene_om /= size;
        branch_om /= size;
        dev_om /= size;

        ofstream nos((name + ".dsomstats").c_str());
        nos << "    \tmean\t\tc-gene\tc-brnch\tc-dev\n";
        nos << fixed << setw(5) << setprecision(2);
        nos << "dS\t" << mean_syn << '\t' << '\t' << gene_syn << '\t' << branch_syn << '\t' << dev_syn << '\n';
        nos << "dN/dS\t" << mean_om << '\t' << '\t' << gene_om << '\t' << branch_om << '\t' << dev_om << '\n';
        cerr << "variance contributions for dS and dN/dS in " << name << ".dsomstats\n";

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
                syn_val[i][j] /= size;
                om_val[i][j] /= size;
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
        devos << "#genename";
        devos << "\tleft_taxon\tright_taxon";
        devos << "\tbranchsynmean\tbranchsynrelvar\tgenesynmean\tgenesynrelvar\tsyn\tsynz";
        devos << "\tbranchommean\tbranchomrelvar\tgeneommean\tgeneomrelvar\tom\tomz";
        devos << "\tdscount\tdsbeta\tempds";
        devos << "\tdncount\tdnbeta\tempdn\tempdnds";
        devos << '\n';
        for (int i=0; i<Ngene; i++) {
            for (int j=0; j<Nbranch; j++)   {
                devos << GetModel()->GetGeneName(i); 
                devos << '\t' << GetModel()->GetTree()->GetLeftTaxon(j);
                devos << '\t' << GetModel()->GetTree()->GetRightTaxon(j);
                devos << '\t' << mean_branchsyn_array[j] << '\t' << branch_syn_relvar[j]; 
                devos << '\t' << mean_genesyn_array[i] << '\t' << gene_syn_relvar[i]; 
                devos << '\t' << syn_val[i][j] << '\t' << syn_z[i][j];
                devos << '\t' << mean_branchom_array[j] << '\t' << branch_om_relvar[j]; 
                devos << '\t' << mean_geneom_array[i] << '\t' << gene_om_relvar[i];
                devos << '\t' << om_val[i][j] << '\t' << om_z[i][j];

                devos << '\t' << GetModel()->GetSynCount(i,j);
                devos << '\t' << GetModel()->GetSynBeta(i,j);
                devos << '\t' << (GetModel()->GetSynCount(i,j)+1) / (GetModel()->GetSynBeta(i,j)+1);
                devos << '\t' << GetModel()->GetNonSynCount(i,j);
                devos << '\t' << GetModel()->GetNonSynBeta(i,j);
                devos << '\t' << (GetModel()->GetNonSynCount(i,j)+1) / (GetModel()->GetNonSynBeta(i,j)+1);
                devos << '\t' << 
                    (GetModel()->GetNonSynCount(i,j)+1) / (GetModel()->GetNonSynBeta(i,j)+1) /
                    (GetModel()->GetSynCount(i,j)+1) * (GetModel()->GetSynBeta(i,j)+1);
                devos << '\n';
            }
        }

        ofstream gos((name + ".postmean.genedsom.tab").c_str());
        gos << "#genename\tgenesynmean\tgenesynrelvar\tgeneommean\tgeneomrelvar\n";
        for (int i=0; i<Ngene; i++) {
            gos << GetModel()->GetGeneName(i) << '\t' << mean_genesyn_array[i] << '\t' << gene_syn_relvar[i] << '\t' << mean_geneom_array[i] << '\t' << gene_om_relvar[i] << '\n';
        }
        cerr << "post mean gene dS and dN/dS in " << name << ".postmean.genedsom.tab\n";

        ofstream bsos((name + ".postmean.branchds.tab").c_str());
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
};

int main(int argc, char *argv[]) {
    int burnin = 0;
    int every = 1;
    int until = -1;
    int ppdev = 0;
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
            } else if (s == "-mixdev")   {
                ppdev = 1;
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
    if (ppdev) {
        sample->ReadMixDev(z, cutoff);
    }
    else    {
        sample->Read();
    }
}
