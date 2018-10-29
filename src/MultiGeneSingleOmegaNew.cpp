#include <cmath>
#include <fstream>
#include "InferenceAppArgParse.hpp"
#include "MultiGeneChain.hpp"
#include "MultiGeneSingleOmegaModelNew.hpp"
#include "SlaveChainDriver.hpp"
#include "components/ChainDriver.hpp"
#include "components/ConsoleLogger.hpp"
#include "components/StandardTracer.hpp"

using namespace std;
using std::unique_ptr;

class MultiGeneSingleOmegaArgParse : public BaseArgParse {
  public:
    MultiGeneSingleOmegaArgParse(ChainCmdLine& cmd) : BaseArgParse(cmd) {}
    ValueArg<string> omega{"", "omega", "omega hyperparameter, expected syntax: <float>,<float>",
        false, "uninf", "string", cmd};
    ValueArg<string> nucrates{"", "nucrates",
        "Possible values are, shared, shrunken, independent(ind)", false, "shrunken", "string",
        cmd};
    ValueArg<string> bl{"", "blmode", "Possible values are, shared, shrunken, independent(ind)",
        false, "shrunken", "string", cmd};

    param_mode_t blmode() { return decode_mode("bl", bl.getValue()); }
    param_mode_t nucmode() { return decode_mode("nucrates", bl.getValue()); }
    omega_param_t omega_param() {
        omega_param_t omega_param;
        if (omega.getValue() != "uninf") {
            std::stringstream ss(omega.getValue());
            std::string item;
            std::getline(ss, item, ',');
            omega_param.hypermean = atof(item.c_str());
            std::getline(ss, item, ',');
            omega_param.hyperinvshape = atof(item.c_str());
            omega_param.variable = false;
        } else {
            omega_param.hypermean = 1.0;
            omega_param.hyperinvshape = 1.0;
            omega_param.variable = true;
        }
        return omega_param;
    }

  private:
    param_mode_t decode_mode(std::string option_name, std::string mode) {
        if (mode == "shared") {
            return shared;
        } else if (mode == "shrunken") {
            return shrunken;
        } else if ((mode == "ind") || (mode == "independent")) {
            return independent;
        } else {
            cerr << "error: does not recognize command after -" << option_name << "\n";
            exit(1);
        }
    }
};

template <class D, class M>
struct AppData {
    unique_ptr<D> chain_driver;
    unique_ptr<M> model;
};

template <class D, class M>
AppData<D, M> load_appdata(ChainCmdLine& cmd, int myid, int nprocs) {
    if (cmd.resume_from_checkpoint()) {
        AppData<D, M> d;
        std::ifstream is = cmd.checkpoint_file();
        d.chain_driver = unique_ptr<D>(new D(is));
        // model = new SingleOmegaModel(is);
        return d;
    } else {
        AppData<D, M> d;
        InferenceAppArgParse app(cmd);
        MultiGeneSingleOmegaArgParse args(cmd);
        cmd.parse();
        d.chain_driver =
            unique_ptr<D>(new D(cmd.chain_name(), app.every.getValue(), app.until.getValue()));
        d.model =
            unique_ptr<M>(new M(app.alignment.getValue(), app.treefile.getValue(), myid, nprocs,
                                args.blmode(), args.nucmode(), args.omega_param()));
        d.model->Allocate();
        d.model->Update();
        return d;
    }
}

int main(int argc, char* argv[]) {
    int myid = 0;
    int nprocs = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


    ChainCmdLine cmd{argc, argv, "MultiGeneSingleOmega", ' ', "0.1"};

    if (!myid) {
        auto d = load_appdata<ChainDriver, MultiGeneSingleOmegaModelMaster>(cmd, myid, nprocs);
        ConsoleLogger console_logger;
        StandardTracer trace(*d.model, cmd.chain_name());
        d.chain_driver->add(*d.model);
        d.chain_driver->add(console_logger);
        d.chain_driver->add(trace);
        d.chain_driver->go();
    } else {
        auto d = load_appdata<SlaveChainDriver, MultiGeneSingleOmegaModelSlave>(cmd, myid, nprocs);
        d.chain_driver->add(*d.model);
        d.chain_driver->go();
    }
    MPI_Finalize();
}
