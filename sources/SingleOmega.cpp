#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "SingleOmegaModel.hpp"
using namespace std;

class SingleOmegaLifecycle : public Lifecycle, public tc::Component {
    SingleOmegaModel* model;

  public:
    SingleOmegaLifecycle() { port("model", &SingleOmegaLifecycle::model); }

    void Init() override {
        model->DeclareModel();  // instead of Allocate
        model->Unfold();
        // cerr << "-- Reset" << endl;
        // Reset(bool force);
        // cerr << "-- initial ln prob = " << model->GetLogProb() << "\n";
        // model->Trace(cerr);
    }

    void EndMove() override {}

    void End() override {}

    // FIXME FIXME FIXME FIXME chain restarting no longer supported!
    // void Open() override {
    //     ifstream is((name + ".param").c_str());
    //     if (!is) {
    //         cerr << "-- Error : cannot find file : " << name << ".param\n";
    //         exit(1);
    //     }
    //     is >> modeltype;
    //     is >> datafile >> treefile;
    //     int tmp;
    //     is >> tmp;
    //     if (tmp) {
    //         cerr << "-- Error when reading model\n";
    //         exit(1);
    //     }
    //     is >> every >> until >> size;

    //     if (modeltype == "SINGLEOMEGA") {
    //         model = new SingleOmegaModel(datafile, treefile);
    //     } else {
    //         cerr << "-- Error when opening file " << name << " : does not recognise model type : " << modeltype << '\n';
    //         exit(1);
    //     }
    //     // model->Allocate();
    //     model->FromStream(is);
    //     model->Update();
    //     model->Unfold();
    //     cerr << size << " points saved, current ln prob = " << model->GetLogProb() << "\n";
    //     model->Trace(cerr);
    // }

    // void Save() override {
    //     ofstream param_os((name + ".param").c_str());
    //     param_os << GetModelType() << '\n';
    //     param_os << datafile << '\t' << treefile << '\n';
    //     param_os << 0 << '\n';
    //     param_os << every << '\t' << until << '\t' << size << '\n';
    //     model->ToStream(param_os);
    // }
};

int main(int argc, char* argv[]) {
    Model model;
    string name;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        cerr << "ERROR: restarting chain not implemented in component version!\n";
        exit(1);
        // SingleOmegaLifecycle* chain = new SingleOmegaLifecycle(name);
        // cerr << "chain " << name << " started\n";
        // chain->start();
        // cerr << "chain " << name << " stopped\n";
        // cerr << chain->GetSize() << " points saved, current ln prob = " << chain->model->GetLogProb() << "\n";
        // chain->model->Trace(cerr);
    }

    // new chain
    else {
        string datafile = "";
        string treefile = "";
        int every = 1;
        int until = -1;

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
                    // } else if (s == "-f") {
                    //     force = 1;
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

        using SOM = SingleOmegaModel;
        using SOMTrace = TraceFile<SOM>;

        model.component<SOM>("model", datafile, treefile);  // remove files form this constructor

        model.component<SingleOmegaLifecycle>("lifecycle").connect<Use<SOM>>("model", "model");

        model.component<ChainDriver>("chaindriver", every, until)
            .connect<Use<ProbModel>>("model", "model")
            .connect<Use<Lifecycle>>("lifecycle", "lifecycle")
            .connect<Use<AbstractTraceFile>>("chainfile", "chainfile")
            .connect<Use<AbstractTraceFile>>("monitorfile", "monitorfile")
            // .connect<Use<AbstractTraceFile>>("paramfile", "paramfile")
            .connect<Use<AbstractTraceFile>>("tracefile", "tracefile")
            .connect<Use<RunToggle>>("runtoggle", "runtoggle");

        model.component<SOMTrace>("chainfile", name + ".chain", &SOM::ToStream, &SOM::ChainHeader)
            .connect<Use<SOM>>("target", "model");

        model.component<SOMTrace>("tracefile", name + ".trace", &SOM::Trace, &SOM::TraceHeader)
            .connect<Use<SOM>>("target", "model");

        model.component<SOMTrace>("monitorfile", name + ".monitor", &SOM::Monitor, nullptr)
            .connect<Use<SOM>>("target", "model");

        // model.component<SOMTrace>>("paramfile", name + ".param");
        model.component<RunToggle>("runtoggle", name);
    }
    Assembly assembly(model);  // instantiating assembly!

    // SingleOmegaLifecycle* chain = new SingleOmegaLifecycle(datafile, treefile, every, until, name, force);
    cerr << "-- Chain " << name << " starting...\n";
    // chain->start();
    assembly.call("chaindriver", "go");
    cerr << "-- Chain " << name << " stopped\n";
    // cerr << chain->GetSize() << "-- Points saved, current ln prob = " << chain->model->GetLogProb() << "\n";
    // chain->model->Trace(cerr);
}
