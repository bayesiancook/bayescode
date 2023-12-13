
#include <cmath>
#include <fstream>
#include "Sample.hpp"
#include "FastCoevolModel.hpp"
#include "DistBranchNodeArray.hpp"
#include "MeanCovMatrix.hpp"
using namespace std;

/**
 * \brief An MCMC sample for FastCoevolModel
 *
 * implements a simple read function, returning the MCMC estimate of the
 * posterior mean and standard deviation of omega=dN/dS
 */

class FastCoevolSample : public Sample {
  private:
    string modeltype;
    string contdatafile, treefile, rootfile;
    string dsomsuffstatfile;
    int wndsmode, wnommode;

  public:
    string GetModelType() override { return modeltype; }

    FastCoevolModel *GetModel() override { return (FastCoevolModel *)model; }

    //! \brief Constructor (file name, burn-in, thinning and upper limit, see
    //! Sample)
    FastCoevolSample(string filename, int inburnin, int inevery, int inuntil)
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
        is >> contdatafile >> treefile >> rootfile;
        is >> dsomsuffstatfile;
        is >> wndsmode >> wnommode;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> chainevery >> chainuntil >> chainsize;

        // make a new model depending on the type obtained from the file
        if (modeltype == "FASTCOEVOLDNDS") {
            model = new FastCoevolModel(contdatafile, treefile, rootfile, dsomsuffstatfile, wndsmode, wnommode);
        } else {
            cerr << "error when opening file " << name << '\n';
            exit(1);
        }

        GetModel()->Allocate();

        // read model (i.e. chain's last point) from <name>.param
        model->FromStream(is);

        // open <name>.chain, and prepare stream and stream iterator
        OpenChainFile();
        // now, size is defined (it is the total number of points with which this
        // Sample object will make all its various posterior averages) all these
        // points can be accessed to (only once) by repeated calls to GetNextPoint()
    }

    void InitializeContrasts(const Link* from, 
            map<const Link*, vector<double>>& contrasts,
            map<const Link*, double>& lengths) {
        if (! from->isLeaf())   {
            contrasts[from] = vector<double>(GetModel()->GetCovMatrix().GetDim(), 0);
            lengths[from] = 0;
            for (const Link* link=from->Next(); link!=from; link=link->Next())  {
                InitializeContrasts(link->Out(), contrasts, lengths);
            }
        }
    }

    void NormalizeContrasts(const Link* from,
            map<const Link*, vector<double>>& contrasts,
            map<const Link*, double>& lengths,
            int size)   {
        if (! from->isLeaf())   {
            lengths[from] /= size;
            vector<double>& cons = contrasts[from];
            for (size_t i=0; i<cons.size(); i++)    {
                cons[i] /= size;
            }
            for (const Link* link=from->Next(); link!=from; link=link->Next())  {
                NormalizeContrasts(link->Out(), contrasts, lengths, size);
            }
        }
    }

    void OutputContrasts(const Link* from,
            map<const Link*, vector<double>>& contrasts,
            map<const Link*, double>& lengths,
            ostream& os)    {
        if (! from->isLeaf())   {
            os << lengths[from];
            vector<double>& cons = contrasts[from];
            for (size_t i=0; i<cons.size(); i++)    {
                os << '\t' << cons[i];
            }
            os << '\n';
            for (const Link* link=from->Next(); link!=from; link=link->Next())  {
                OutputContrasts(link->Out(), contrasts, lengths, os);
            }
        }
    }

    void ReadIndependentContrasts() {
        cerr << size << " points to read\n";

        map<const Link*, vector<double>> contrasts;
        map<const Link*, double> lengths;
        InitializeContrasts(GetModel()->GetTree().GetRoot(), contrasts, lengths);

        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();
            GetModel()->GetProcess().GetIndependentContrasts(contrasts, lengths);
        }
        cerr << '\n';

        NormalizeContrasts(GetModel()->GetTree().GetRoot(), contrasts, lengths, size);
        ofstream os((name + ".ic").c_str());
        OutputContrasts(GetModel()->GetTree().GetRoot(), contrasts, lengths, os);
        cerr << "independent contrasts tabulated in " << name << ".ic\n";
    }

    //! \brief computes the posterior mean estimate (and the posterior standard
    //! deviation) of omega
    void Read() {
        cerr << size << " points to read\n";

        DistBranchNodeArray meansynrate(GetModel()->GetTree());
        DistBranchNodeArray meanomega(GetModel()->GetTree());
        // DistBranchArray<double> timetree(GetModel()->GetTree());

		int dim = GetModel()->GetCovMatrix().GetDim();
		MeanCovMatrix  mat(dim);

        double pvar_syn = 0;
        double pvar_om = 0;

        for (int i=0; i<size; i++) {
            cerr << '.';
            GetNextPoint();
            GetModel()->Update();
            meansynrate.AddFromChrono(GetModel()->GetChronogram(), GetModel()->GetProcess(), 0);
            meanomega.AddFromChrono(GetModel()->GetChronogram(), GetModel()->GetProcess(), 1);
            // timetree.Add(GetModel()->GetChronogram());
			mat.Add(GetModel()->GetCovMatrix());
            pvar_syn += GetModel()->GetLongTermSynPropVar();
            pvar_om += GetModel()->GetLongTermOmegaPropVar();
        }
        cerr << '\n';

        pvar_syn /= size;
        pvar_om /= size;
        ofstream vos((name + ".propvar").c_str());
        vos << "proportion of variance contributed by long term (Brownian) trends\n";
        vos << "dS\t" << pvar_syn << '\n';
        vos << "dN/dS\t" << pvar_om << '\n';
        cerr << "proportion of variance contributed by long term trends in " << name << ".propvar\n";
    
        meansynrate.Sort();
        ofstream sos((name + ".postmeands.tre").c_str());
        meansynrate.MedianToStream(sos);
        cerr << "postmean dS tree in " << name << ".postmeands.tre\n"; 

        ofstream node_sos((name + ".postmeannodeds.tab").c_str());
        meansynrate.TabulateNodeMedianToStream(node_sos);
        cerr << "tabulated node dS median values in " << name << ".postmeannodeds.tab\n"; 

        ofstream branch_sos((name + ".postmeanbranchds.tab").c_str());
        meansynrate.TabulateBranchMedianToStream(branch_sos);
        cerr << "tabulated branch dS median values in " << name << ".postmeanbranchds.tab\n"; 

        meanomega.Sort();
        ofstream omos((name + ".postmeanomega.tre").c_str());
        meanomega.MedianToStream(omos);
        cerr << "postmean omega tree in " << name << ".postmeanomega.tre\n"; 

        ofstream node_omos((name + ".postmeannodeomega.tab").c_str());
        meanomega.TabulateNodeMedianToStream(node_omos);
        cerr << "tabulated node omega median values in " << name << ".postmeannodeomega.tab\n"; 

        ofstream branch_omos((name + ".postmeanbranchomega.tab").c_str());
        meanomega.TabulateBranchMedianToStream(branch_omos);
        cerr << "tabulated branch omega median values in " << name << ".postmeanbranchomega.tab\n"; 

        // timetree.Sort();
        ofstream tos((name + ".postmeanchrono.tre").c_str());
        meanomega.MedianToStream(tos, false);
        cerr << "postmean timetree in " << name << ".postmeanchrono.tre\n"; 

		mat.Normalize();
		ofstream mos((name + ".cov").c_str());
		mos << "entries are in the following order:\n";
		GetModel()->PrintEntries(mos);
		mos << '\n';
		// mat.SetLatex(tex);
		mos << mat;
		cerr << "covariance matrix in " << name << ".cov\n";
		cerr << '\n';
    }
};

int main(int argc, char *argv[]) {
    int burnin = 0;
    int every = 1;
    int until = -1;
    int ppred = 0;
    int ic = 0;

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
            } else if (s == "-ic")  {
                ic = 1;
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

    FastCoevolSample *sample = new FastCoevolSample(name, burnin, every, until);
    if (ppred) {
        sample->PostPred();
    } else if (ic)  {
        sample->ReadIndependentContrasts();
    } else {
        sample->Read();
    }
}
