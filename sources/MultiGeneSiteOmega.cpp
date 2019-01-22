#include <cmath>
#include <fstream>
#include "MultiGeneChain.hpp"
#include "MultiGeneSiteOmegaModel.hpp"

using namespace std;

MPI_Datatype Propagate_arg;

/**
 * \brief Chain object for running an MCMC under MultiGeneSiteOmegaModel
 */

class MultiGeneSiteOmegaChain : public MultiGeneChain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int writegenedata;
    int blmode, nucmode, omegamode;
    double omegameanhypermean, omegameanhyperinvshape;
    double omegainvshapehypermean, omegainvshapehyperinvshape;

  public:
    MultiGeneSiteOmegaModel *GetModel() {
        return static_cast<MultiGeneSiteOmegaModel *>(model);
    }

    string GetModelType() override { return modeltype; }

    //! \brief constructor for a new MCMC
    //!
    //! \param indatafile: name of file contanining sequence alignment
    //! \param intreefile: name of file contaning tree
    //! \param inevery: thinning factor
    //! \param inuntil: maximum MCMC sample size
    //! \param name: base name for all files related to this MCMC run
    //! \param force: overwrite existing files with same name
    //! \param inmyid, int innprocs: process id and total number of MPI processes
    MultiGeneSiteOmegaChain(string indatafile, string intreefile, int inblmode, int innucmode, int inomegamode,
                              double inomegameanhypermean, double inomegameanhyperinvshape, 
                              double inomegainvshapehypermean, double inomegainvshapehyperinvshape, 
                              int inevery, int inuntil, int inwritegenedata,
                              string inname, int force, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs),
          modeltype("MULTIGENESITEOMEGA"),
          datafile(indatafile),
          treefile(intreefile) {
        blmode = inblmode;
        nucmode = innucmode;
        omegamode = inomegamode;
        omegameanhypermean = inomegameanhypermean;
        omegameanhyperinvshape = inomegameanhyperinvshape;
        omegainvshapehypermean = inomegainvshapehypermean;
        omegainvshapehyperinvshape = inomegainvshapehyperinvshape;
        every = inevery;
        until = inuntil;
        writegenedata = inwritegenedata;
        name = inname;
        New(force);
    }

    //! \brief constructor for opening and restarting an already existing chain
    MultiGeneSiteOmegaChain(string filename, int inmyid, int innprocs)
        : MultiGeneChain(inmyid, innprocs) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new MultiGeneSiteOmegaModel(datafile, treefile, myid, nprocs);
        GetModel()->SetAcrossGenesModes(blmode,nucmode,omegamode);
        GetModel()->SetOmegaHyperParameters(omegameanhypermean, omegameanhyperinvshape, omegainvshapehypermean, omegainvshapehyperinvshape);
        if (!myid) {
            cerr << "allocate\n";
        }
        GetModel()->Allocate();
        if (!myid) {
            cerr << "update\n";
        }
        GetModel()->Update();
        Reset(force);
        if (!myid) {
            model->Trace(cerr);
        }
    }

    void Open() override {
        ifstream is((name + ".param").c_str());
        if (!is) {
            cerr << "error : cannot find file : " << name << ".param\n";
            exit(1);
        }
        is >> modeltype;
        is >> datafile >> treefile;
        is >> writegenedata;
        is >> blmode >> nucmode >> omegamode;
        is >> omegameanhypermean >> omegameanhyperinvshape;
        is >> omegainvshapehypermean >> omegainvshapehyperinvshape;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "MULTIGENESITEOMEGA") {
            model = new MultiGeneSiteOmegaModel(datafile, treefile, myid, nprocs);
            GetModel()->SetAcrossGenesModes(blmode,nucmode,omegamode);
            GetModel()->SetOmegaHyperParameters(omegameanhypermean, omegameanhyperinvshape, omegainvshapehypermean, omegainvshapehyperinvshape);
        } else {
            cerr << "error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        if (!myid) {
            cerr << "allocate\n";
        }
        GetModel()->Allocate();
        model->FromStream(is);
        if (!myid) {
            cerr << "update\n";
        }
        model->Update();
        if (!myid) {
            cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
            model->Trace(cerr);
        }
    }

    void Save() override {
        if (!myid) {
            ofstream param_os((name + ".param").c_str());
            param_os << GetModelType() << '\n';
            param_os << datafile << '\t' << treefile << '\n';
            param_os << writegenedata << '\n';
            param_os << blmode << '\t' << nucmode << '\t' << omegamode << '\n';
            param_os << omegameanhypermean << '\t' << omegameanhyperinvshape << '\n';
            param_os << omegainvshapehypermean << '\t' << omegainvshapehyperinvshape << '\n';
            param_os << 0 << '\n';
            param_os << every << '\t' << until << '\t' << size << '\n';
            GetModel()->MasterToStream(param_os);
        } else {
            GetModel()->SlaveToStream();
        }
    }

    void MakeFiles(int force) override {
        if (myid) {
            cerr << "error: in specialized makefiles\n";
            exit(1);
        }
        MultiGeneChain::MakeFiles(force);
        if (writegenedata)  {
            ofstream os((name + ".geneom").c_str());
            if (writegenedata == 2) {
                ofstream os((name + ".siteom").c_str());
            }
        }
    }

    void SavePoint() override {
        MultiGeneChain::SavePoint();
        if (writegenedata == 1)  {
            if (!myid) {
                ofstream os((name + ".geneom").c_str(), ios_base::app);
                GetModel()->TraceOmega(os);
            }
        }
        if (writegenedata == 2) {
            if (!myid) {
                ofstream os((name + ".siteom").c_str(), ios_base::app);
                GetModel()->MasterTraceSiteOmega(os);
            } else {
                GetModel()->SlaveTraceSiteOmega();
            }
        }
    }
};

int main(int argc, char *argv[]) {
    int myid = 0;
    int nprocs = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    MultiGeneSiteOmegaChain *chain = 0;
    string name = "";

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new MultiGeneSiteOmegaChain(name, myid, nprocs);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        int force = 1;
        int every = 1;
        int until = -1;
        int blmode = 1;
        int nucmode = 1;
        int omegamode = 1;
        double omegameanhypermean = 1.0;
        double omegameanhyperinvshape = 1.0;
        double omegainvshapehypermean = 1.0;
        double omegainvshapehyperinvshape = 1.0;
        int writegenedata = 1;

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
                } else if (s == "-f") {
                    force = 1;
                } else if (s == "-omega") {
                    omegamode = 0;
                    i++;
                    string tmp = argv[i];
                    if (tmp != "uninf") {
                        omegameanhypermean = atof(argv[i]);
                        i++;
                        omegameanhyperinvshape = atof(argv[i]);
                        i++;
                        omegainvshapehypermean = atof(argv[i]);
                        i++;
                        omegainvshapehyperinvshape = atof(argv[i]);
                    }
                } else if (s == "-nucrates") {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "shared") {
                        nucmode = 2;
                    } else if (tmp == "shrunken") {
                        nucmode = 1;
                    } else if ((tmp == "ind") || (tmp == "independent")) {
                        nucmode = 0;
                    } else {
                        cerr << "error: does not recongnize command after -nucrates\n";
                        exit(1);
                    }
                } else if (s == "-bl") {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "shared") {
                        blmode = 2;
                    } else if (tmp == "shrunken") {
                        blmode = 1;
                    } else if ((tmp == "ind") || (tmp == "independent")) {
                        blmode = 0;
                    } else {
                        cerr << "error: does not recongnize command after -bl\n";
                        exit(1);
                    }
                } else if (s == "-g") {
                    writegenedata = 0;
                } else if (s == "+g") {
                    writegenedata = 1;
                } else if (s == "+G") {
                    writegenedata = 2;
                } else if ((s == "-x") || (s == "-extract")) {
                    i++;
                    if (i == argc) throw(0);
                    every = atoi(argv[i]);
                    i++;
                    if (i == argc) throw(0);
                    until = atoi(argv[i]);
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
            cerr << "globom -d <alignment> -t <tree> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new MultiGeneSiteOmegaChain(datafile, treefile, blmode, nucmode, omegamode,
                                              omegameanhypermean, omegameanhyperinvshape, omegainvshapehypermean, omegainvshapehyperinvshape,
                                              every, until, writegenedata, name, force, myid,
                                              nprocs);
    }

    if (!myid) {
        cerr << "chain " << name << " started\n";
    }
    chain->Start();
    if (!myid) {
        cerr << "chain " << name << " stopped\n";
        cerr << chain->GetSize()
             << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        chain->GetModel()->Trace(cerr);
    }

    MPI_Finalize();
}
