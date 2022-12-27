#include <cmath>
#include <fstream>
#include "Sample.hpp"
#include "SingleOmegaModel.hpp"
using namespace std;

/**
 * \brief An MCMC sample for SingleOmegaModel
 *
 * implements a simple read function, returning the MCMC estimate of the
 * posterior mean and standard deviation of omega=dN/dS
 */

class SingleOmegaSample : public Sample {
  private:
    string modeltype;
    string datafile;
    string treefile;
    int fixbl;

  public:
    string GetModelType() override { return modeltype; }

    SingleOmegaModel *GetModel() override { return (SingleOmegaModel *)model; }

    //! \brief Constructor (file name, burn-in, thinning and upper limit, see
    //! Sample)
    SingleOmegaSample(string filename, int inburnin, int inevery, int inuntil)
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
        fixbl = 0;
        int check;
        is >> check;
        if (check) {
            is >> fixbl;
            is >> check;
            if (check)  {
                cerr << "-- Error when reading model\n";
                exit(1);
            }
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "SINGLEOMEGA") {
            model = new SingleOmegaModel(datafile, treefile);
        } else {
            cerr << "error when opening file " << name << '\n';
            exit(1);
        }

        GetModel()->SetFixBL(fixbl);
        GetModel()->Allocate();

        // read model (i.e. chain's last point) from <name>.param
        model->FromStream(is);

        // open <name>.chain, and prepare stream and stream iterator
        OpenChainFile();
        // now, size is defined (it is the total number of points with which this
        // Sample object will make all its various posterior averages) all these
        // points can be accessed to (only once) by repeated calls to GetNextPoint()
    }

    void ReadSiteSuffStat() {
        cerr << size << " points to read\n";

        dSOmegaPathSuffStatArray array(GetModel()->GetNsite());
        array.Clear();

        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();
            GetModel()->AdddSOmegaPathSuffStat(array);
        }
        cerr << '\n';
        array.Normalize(1.0/size);
        ofstream os((name + ".sitedsomss").c_str());
        os << "site\tMs\tMn\tLs\tLn\tom\n";
        for (int i=0; i<array.GetSize(); i++)   {
            os << i << '\t' << array[i] << '\n';
            // os << i << '\t' << array[i] << '\t' << array[i].GetdNdS(0.1) << '\n';
        }
        cout << "site dsom suffstat written in : " << name << ".sitedsomss\n";
    }

    void ReadNodePathSuffStat() {
        cerr << size << " points to read\n";
        RelativePathSuffStatNodeArray array(*GetModel()->GetTree(), GetModel()->GetCodonStateSpace()->GetNstate());
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();
            GetModel()->AddNodePathSuffStat(array);
        }
        cerr << '\n';
        array.Normalize(1.0/size);
        ofstream os((name + ".meannodepathsuffstat").c_str());
        os << array << '\n';
        cerr << "node path suffstats in " << name << ".meannodepathsuffstat\n";
    }

    void ReadBranchSuffStat()   {
        cerr << size << " points to read\n";
        dSOmegaPathSuffStatBranchArray array(*GetModel()->GetTree());
        GCConsdSOmegaPathSuffStatBranchArray gcconsarray(*GetModel()->GetTree());
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();
            GetModel()->AdddSOmegaPathSuffStat(array);
            GetModel()->AddGCConsdSOmegaPathSuffStat(gcconsarray);
        }
        cerr << '\n';
        array.Normalize(1.0/size);
        gcconsarray.Normalize(1.0/size);

        /*
        ofstream gos((name + ".meanbranchdsomsuffstat").c_str());
        gos << array << '\n';
        cerr << "global dsom path suffstats in " << name << ".meanbranchdsomsuffstat\n";
        ofstream tgos((name + ".meanbranchdsomsuffstat.tre").c_str());
        array.BranchToNewick(tgos);
        cerr << "newick format: global dsom path suffstats in " << name << ".meanbranchdsomsuffstat.tre\n";

        ofstream gcgos((name + ".meanbranchgcconsdsomsuffstat").c_str());
        gcgos << gcconsarray << '\n';
        cerr << "global gc-cons dsom path suffstats in " << name << ".meanbranchgcconsdsomsuffstat\n";
        ofstream gctgos((name + ".meanbranchgcconsdsomsuffstat.tre").c_str());
        // gcconsarray.BranchToNewick(gctgos);
        // cerr << "newick format: global gc-cons dsom path suffstats in " << name << ".meanbranchgcconsdsomsuffstat.tre\n";
        */

        ofstream syn_count_os((name + ".counts_dS.dnd").c_str());
        array.BranchToNewickSynCount(syn_count_os);
        ofstream syn_beta_os((name + ".counts_dS_norm.dnd").c_str());
        array.BranchToNewickSynBeta(syn_beta_os);
        ofstream nonsyn_count_os((name + ".counts_dN.dnd").c_str());
        array.BranchToNewickNonSynCount(nonsyn_count_os);
        ofstream nonsyn_beta_os((name + ".counts_dN_norm.dnd").c_str());
        array.BranchToNewickNonSynBeta(nonsyn_beta_os);
        cerr << "dsom path suffstats in " << name << ".counts_dS.dnd, counts_dS_norm.dnd, counts_dN.dnd, counts_dN_norm.dnd\n";

        ofstream gccons_syn_count_os((name + ".gccons.counts_dS.dnd").c_str());
        gcconsarray.BranchToNewickSynCount(gccons_syn_count_os);
        ofstream gccons_syn_beta_os((name + ".gccons.counts_dS_norm.dnd").c_str());
        gcconsarray.BranchToNewickSynBeta(gccons_syn_beta_os);
        ofstream gccons_nonsyn_count_os((name + ".gccons.counts_dN.dnd").c_str());
        gcconsarray.BranchToNewickNonSynCount(gccons_nonsyn_count_os);
        ofstream gccons_nonsyn_beta_os((name + ".gccons.counts_dN_norm.dnd").c_str());
        gcconsarray.BranchToNewickNonSynBeta(gccons_nonsyn_beta_os);
        cerr << "dsom path suffstats in " << name << ".counts_dS.dnd, counts_dS_norm.dnd, counts_dN.dnd, counts_dN_norm.dnd\n";
    }

    //! \brief computes the posterior mean estimate (and the posterior standard
    //! deviation) of omega
    void Read() {
        cerr << size << " points to read\n";

        double meanomega = 0;
        double varomega = 0;

        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            double om = GetModel()->GetOmega();
            meanomega += om;
            varomega += om * om;
        }
        cerr << '\n';
        meanomega /= size;
        varomega /= size;
        varomega -= meanomega * meanomega;

        cout << "posterior mean omega : " << meanomega << '\t' << sqrt(varomega) << '\n';
    }
};

int main(int argc, char *argv[]) {
    int burnin = 0;
    int every = 1;
    int until = -1;
    int ppred = 0;
    int sitess = 0;
    int branchss = 0;
    int nodepathss = 0;

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
            } else if (s == "-sitedsomss")    {
                sitess = 1;
            } else if (s == "-branchdsomss")    {
                branchss = 1;
            } else if (s == "-nodepathss")  {
                nodepathss = 1;
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

    SingleOmegaSample *sample = new SingleOmegaSample(name, burnin, every, until);
    if (ppred) {
        sample->PostPred();
    } else if (sitess)  {
        sample->ReadSiteSuffStat();
    } else if (branchss)  {
        sample->ReadBranchSuffStat();
    } else if (nodepathss)  {
        sample->ReadNodePathSuffStat();
    } else {
        sample->Read();
    }
}
