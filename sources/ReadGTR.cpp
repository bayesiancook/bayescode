#include <cmath>
#include <fstream>
#include "Sample.hpp"
#include "GTRModel.hpp"
using namespace std;

/**
 * \brief An MCMC sample for GTRModel
 *
 * implements a simple read function, returning the MCMC estimate of the
 * posterior mean and standard deviation of omega=dN/dS
 */

class GTRSample : public Sample {
  private:
    string modeltype;
    string datafile;
    string treefile;
    int fixbl;
    int fixrr;

  public:
    string GetModelType() override { return modeltype; }

    GTRModel *GetModel() override { return (GTRModel *)model; }

    //! \brief Constructor (file name, burn-in, thinning and upper limit, see
    //! Sample)
    GTRSample(string filename, int inburnin, int inevery, int inuntil)
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
        is >> fixbl >> fixrr;
        int check;
        is >> check;
        if (check) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "GTR") {
            model = new GTRModel(datafile, treefile);
        } else {
            cerr << "error when opening file " << name << '\n';
            exit(1);
        }

        GetModel()->SetFixBL(fixbl);
        GetModel()->Allocate();
        if (fixrr)  {
            GetModel()->SetLG();
        }

        // read model (i.e. chain's last point) from <name>.param
        model->FromStream(is);

        // open <name>.chain, and prepare stream and stream iterator
        OpenChainFile();
        // now, size is defined (it is the total number of points with which this
        // Sample object will make all its various posterior averages) all these
        // points can be accessed to (only once) by repeated calls to GetNextPoint()
    }

    void ReadNodePathSuffStat() {
        cerr << size << " points to read\n";
        RelativePathSuffStatNodeArray array(*GetModel()->GetTree(), GetModel()->GetStateSpace()->GetNstate());
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

    void ReadSiteNodePathSuffStat() {
        cerr << size << " points to read\n";
        PathSuffStatBidimArray array(GetModel()->GetTree()->GetNbranch(), GetModel()->GetNsite(), GetModel()->GetStateSpace()->GetNstate());
        BranchAllocationSystem alloc(*GetModel()->GetTree(), GetModel()->GetTree()->GetNbranch());
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();
            GetModel()->AddPathSuffStat(array, alloc);
        }
        cerr << '\n';
        array.Normalize(1.0/size);

        ofstream os((name + ".suffstat").c_str());
        const double aa_deg[] = {4.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0, 3.0, 6.0, 6.0, 1.0, 2.0, 4.0, 2.0, 4.0, 4.0, 4.0, 4.0, 1.0, 2.0};
        int nstate = GetModel()->GetStateSpace()->GetNstate();
        int nsite = GetModel()->GetNsite();
        int nbranch = GetModel()->GetTree()->GetNbranch();

        // header
        os << nsite << '\t' << nbranch << '\n';

        // root degeneracy (1 line)
        for (int j=0; j<nstate; j++)    {
            os << aa_deg[j] << '\t';
        }
        os << '\n';

        // root post probs (nsite lines)
        for (int i=0; i<nsite; i++)    {
            for (int j=0; j<nstate; j++)    {
                os << array.GetVal(0,i).GetRootCount(j) << '\t';
            }
            os << '\n';
        }

        // state post probs across branches (nbranch*nsite lines)
        for (int i=0; i<nsite; i++)    {
            for (int b=0; b<nbranch; b++)   {
                double tot = 0;
                for (int j=0; j<nstate; j++)    {
                    tot += array.GetVal(b,i).GetWaitingTime(j);
                }
                if (tot > 0)   {
                    for (int j=0; j<nstate; j++)    {
                        os << array.GetVal(b,i).GetWaitingTime(j) / tot << '\t';
                    }
                }
                else    {
                    for (int j=0; j<nstate; j++)    {
                        os << 0 << '\t';
                    }
                }
                os << '\n';
            }
        }

        // state counts
        for (int i=0; i<nsite; i++)    {
            for (int b=0; b<nbranch; b++)   {
                for (int j=0; j<nstate; j++)    {
                    for (int k=0; k<nstate; k++)    {
                        os << array.GetVal(b,i).GetPairCount(j,k) << '\t';
                    }
                    os << '\n';
                }
            }
        }

        // state mu factors
        for (int i=0; i<nsite; i++)    {
            for (int b=0; b<nbranch; b++)   {
                for (int j=0; j<nstate; j++)    {
                    for (int k=0; k<nstate; k++)    {
                        if (k != j) {
                            os << array.GetVal(b,i).GetWaitingTime(j) * aa_deg[k] << '\t';
                        }
                        else {
                            os << 0 << '\t';
                        }
                    }
                    os << '\n';
                }
            }
        }

        cerr << "path suffstats in " << name << ".suffstat\n";
    }

    //! \brief computes the posterior mean estimate (and the posterior standard
    //! deviation) of omega
    void Read() {
        cerr << size << " points to read\n";
        for (int i = 0; i < size; i++) {
            cerr << '.';
            GetNextPoint();
        }
        cerr << '\n';
    }
};

int main(int argc, char *argv[]) {
    int burnin = 0;
    int every = 1;
    int until = -1;
    int ppred = 0;
    int nodepathss = 0;
    int pathss = 0;

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
            } else if (s == "-nodepathss")  {
                nodepathss = 1;
            } else if (s == "-ss")  {
                pathss = 1;
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

    GTRSample *sample = new GTRSample(name, burnin, every, until);
    if (ppred) {
        sample->PostPred();
    } else if (nodepathss)  {
        sample->ReadNodePathSuffStat();
    } else if (pathss)  {
        sample->ReadSiteNodePathSuffStat();
    } else {
        sample->Read();
    }
}
