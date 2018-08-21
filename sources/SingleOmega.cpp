#include <cmath>
#include <fstream>
#include "BaseArgParse.hpp"
#include "Chain.hpp"
#include "ChainDriver.hpp"
#include "SingleOmegaModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under SingleOmegaModel
 */

class SingleOmegaChain : public Chain {
  private:
    // Chain parameters
    string modeltype;
    string datafile, treefile;

  public:
    //! constructor for a new chain: datafile, treefile, saving frequency, final
    //! chain size, chain name and overwrite flag -- calls New
    SingleOmegaChain(string indatafile, string intreefile, int inevery, int inuntil, string inname,
                     int force)
        : modeltype("SINGLEOMEGA"), datafile(indatafile), treefile(intreefile) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    //! constructor for re-opening an already existing chain from file -- calls
    //! Open
    SingleOmegaChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new SingleOmegaModel(datafile, treefile);
        GetModel()->Allocate();
        GetModel()->Update();
        cerr << "-- Reset" << endl;
        Reset(force);
        cerr << "-- initial ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Open() override {
        ifstream is((name + ".param").c_str());
        if (!is) {
            cerr << "-- Error : cannot find file : " << name << ".param\n";
            exit(1);
        }
        is >> modeltype;
        is >> datafile >> treefile;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "SINGLEOMEGA") {
            model = new SingleOmegaModel(datafile, treefile);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        GetModel()->Allocate();
        model->FromStream(is);
        model->Update();
        cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Save() override {
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datafile << '\t' << treefile << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }

    //! return the model, with its derived type (unlike ProbModel::GetModel)
    SingleOmegaModel *GetModel() { return static_cast<SingleOmegaModel *>(model); }

    //! return model type
    string GetModelType() override { return modeltype; }
};

class SingleOmegaArgParse : public BaseArgParse {
  public:
    SingleOmegaArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd) {}
    ValueArg<string> alignment{"a",      "alignment", "Alignment file (PHYLIP)", true, "",
                               "string", cmd};
    ValueArg<string> treefile{"t", "tree", "Tree file (NHX)", true, "", "string", cmd};
    ValueArg<int> every{"e",   "every", "Number of iterations between two traces", false, 1,
                        "int", cmd};
    ValueArg<int> until{"u",   "until", "Maximum number of (saved) iterations (-1 means unlimited)",
                        false, -1,      "int",
                        cmd};
    SwitchArg force{"f", "force", "Overwrite existing output files", cmd};
};

class ConsoleLogger : public ChainComponent {
  public:
    void start() override { cout << "Started\n"; }
    void move(int i) override { cout << "Move " << i << "\n"; }
    void savepoint(int i) override { cout << "Savepoint " << i << "\n"; }
    void end() override { cout << "Ended\n"; }
};

class ChainCheckpoint : public ChainComponent {
    std::string filename;
    std::function<void(std::ostream &)> serialize_model;
    ChainDriver &cd;

  public:
    template <class T>
    ChainCheckpoint(std::string filename, ChainDriver &cd, T &model)
        : filename(filename),
          serialize_model([&model](std::ostream &os) { model.ToStream(os); }),
          cd(cd) {}

    void savepoint(int) override {
        std::ofstream os{filename};
        cd.serialize(os);
        os << "\n";
        serialize_model(os);
    }
};

template<class T, class... Args>
unique_ptr<T> make_unique(Args&&... args) { return unique_ptr<T>(new T(std::forward<Args>(args)...)); }

int main(int argc, char *argv[]) {
    ChainCmdLine cmd{argc, argv, "SingleOmega", ' ', "0.1"};

    unique_ptr<ChainDriver> chain_driver = nullptr;
    unique_ptr<SingleOmegaModel> model = nullptr;

    if(cmd.resume_from_checkpoint()) {
        std::ifstream is = cmd.checkpoint_file();
        chain_driver = make_unique<ChainDriver>(is);
        model = make_unique<SingleOmegaModel>(is);
    }
    else {
        SingleOmegaArgParse args(cmd);
        cmd.parse();
        chain_driver = make_unique<ChainDriver>(cmd.chain_name(), args.every.getValue(),
                                                args.until.getValue());
        model = make_unique<SingleOmegaModel>(args.alignment.getValue(), args.treefile.getValue());
    }

    ConsoleLogger console_logger;
    ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", *chain_driver, *model);
    chain_driver->add(*model);
    chain_driver->add(console_logger);
    chain_driver->add(chain_checkpoint);
    chain_driver->go();
}
