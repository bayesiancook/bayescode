#include <cmath>
#include <fstream>
#include <sstream>
#include "Chain.hpp"
#include "GTRModel.hpp"
#include "MeanPathSuffStat.hpp"

using namespace std;


/**
 * \brief Chain object for running an MCMC under GTRModel
 */

class GTRChain : public Chain {
  private:
    // Chain parameters
    string modeltype;
    string datafile, treefile;

  public:
    //! constructor for a new chain: datafile, treefile, saving frequency, final
    //! chain size, chain name and overwrite flag -- calls New
    GTRChain(string indatafile, string intreefile, string inname)
        : modeltype("GTR"), datafile(indatafile), treefile(intreefile) {
        name = inname;
        New(1);
    }

    void New(int force) override {
        model = new GTRModel(datafile, treefile);
        GetModel()->SetFixBL(1);
        GetModel()->Allocate();
        GetModel()->SetLG();
        GetModel()->Update();
        GetModel()->TraceHeader(cerr);
        GetModel()->Trace(cerr);
    }

    void Open() override {}
    void Save() override {}

    //! return the model, with its derived type (unlike ProbModel::GetModel)
    GTRModel *GetModel() { return static_cast<GTRModel *>(model); }

    //! return model type
    string GetModelType() override { return modeltype; }

    void ReadSiteNodePathSuffStat(int nrep) {
        MeanPathSuffStatBidimArray array(GetModel()->GetTree()->GetNbranch(), GetModel()->GetNsite(), GetModel()->GetStateSpace()->GetNstate());
        BranchAllocationSystem alloc(*GetModel()->GetTree(), GetModel()->GetTree()->GetNbranch());
        for (int rep=0; rep<nrep; rep++)    {
            cerr << '.';
            GetModel()->Update();
            GetModel()->AddPathSuffStat(array, alloc);
        }
        cerr << '\n';
        array.Normalize(1.0/nrep);

        const double aa_deg[] = {4.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0, 3.0, 6.0, 6.0, 1.0, 2.0, 4.0, 2.0, 4.0, 4.0, 4.0, 4.0, 1.0, 2.0};
        int nsite = GetModel()->GetNsite();
        int nbranch = GetModel()->GetTree()->GetNbranch();

        ofstream os((name + ".suffstat.summary").c_str());
        os << "nsite" << '\t' << nsite << '\n';
        os << "ntaxa" << '\t' << GetModel()->GetNtaxa() << '\n';
        os << "nbranch" << '\t' << nbranch << '\n';
        os << "nstate" << '\t' << 20 << '\n';

        ofstream ppos((name + ".suffstat.branch_state_probs").c_str());
        auto pp_vec = array.GetAllPostProbs();
        ppos.write((char*)pp_vec.data(), pp_vec.size()*sizeof(float));

        ofstream cos((name + ".suffstat.branch_counts").c_str());
        auto c_vec = array.GetAllPairCounts();
        cos.write((char*)c_vec.data(), c_vec.size()*sizeof(float));

        ofstream mos((name + ".suffstat.branch_eff_targets").c_str());
        auto m_vec = array.GetAllWeightedWaitingTimes(aa_deg);
        mos.write((char*)m_vec.data(), m_vec.size()*sizeof(float));

        cerr << "suffstats in " << name << ".suffstat.*\n";
    }

};

int main(int argc, char *argv[]) {
    string name = "";
    GTRChain *chain = 0;

    string datafile = "";
    string treefile = "";
    name = "";
    int nrep = 10;

    try {
        if (argc == 1) {
            throw(0);
        }

        int i = 1;
        while (i < argc) {
            string s = argv[i];

            if (s == "-d") {
                i++;
                datafile = argv[i];
            } else if ((s == "-t") || (s == "-T")) {
                i++;
                treefile = argv[i];
            } else if (s == "-nrep")  {
                i++;
                nrep = atoi(argv[i]);
            } else {
                if (i != (argc - 1)) {
                    throw(0);
                }
                name = argv[i];
            }
            i++;
        }
        if ((datafile == "") || (treefile == "") || (name == "")) {
            throw(0);
        }
    } catch (...) {
        cerr << "fastmapgtr -d <alignment> -t <tree> -nrep <nrep> <name> \n";
        cerr << '\n';
        exit(1);
    }

    chain = new GTRChain(datafile, treefile, name);

    chain->ReadSiteNodePathSuffStat(nrep);
}

