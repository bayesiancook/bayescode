#include <cmath>
#include <fstream>
#include "MultiGeneChain.hpp"
#include "MultiGeneSingleOmegaModelNew.hpp"
#include "InferenceAppArgParse.hpp"
#include "components/ChainDriver.hpp"
#include "SlaveChainDriver.hpp"

using namespace std;
using std::unique_ptr;

class MultiGeneSingleOmegaArgParse : public BaseArgParse {
  public:
    MultiGeneSingleOmegaArgParse(ChainCmdLine &cmd) : BaseArgParse(cmd), app(cmd) {}
    InferenceAppArgParse app;
    ValueArg<string> omega{"","omega","omega hyperparameter, expected syntax: <float>,<float>",false,"uninf","string",cmd};
    ValueArg<string> nucrates{"","nucrates","Possible values are, shared, shrunken, independent(ind)",false,"shrunken","string",cmd};
    ValueArg<string> bl{"","blmode","Possible values are, shared, shrunken, independent(ind)",false,"shrunken","string",cmd};

    int blmode() { return decode_mode("bl", bl.getValue()); }
    int nucmode() { return decode_mode("nucrates", bl.getValue()); }
    int omegamode(double& omegahypermean, double& omegahyperinvshape) {
        if (omega.getValue() != "uninf") {
            std::stringstream ss(omega.getValue());
            std::string item;
            std::getline(ss, item, ',');
            omegahypermean = atof(item.c_str());
            std::getline(ss, item, ',');
            omegahyperinvshape = atof(item.c_str());
            return 0;
        }
        else {
            omegahypermean = 1.0;
            omegahyperinvshape = 1.0;
            return 1;
        }
    }

private:
    int decode_mode(std::string option_name, std::string mode) {
        if (mode == "shared") {
            return 2;
        } else if (mode == "shrunken") {
            return 1;
        } else if ((mode == "ind") || (mode == "independent")) {
            return 0;
        } else {
            cerr << "error: does not recongnize command after -" << option_name << "\n";
            exit(1);
        }
    }
};

template<class D, class M>
struct AppData {
    unique_ptr<D> chain_driver;
    unique_ptr<M> model;
};

template<class D, class M>
AppData<D, M> load_appdata(ChainCmdLine& cmd, int myid, int nprocs) {

    if (cmd.resume_from_checkpoint()) {
        AppData<D, M> d;
        std::ifstream is = cmd.checkpoint_file();
        d.chain_driver = unique_ptr<D>(new D(is));
        // model = new SingleOmegaModel(is);
        return d;
    } else {
        AppData<D, M> d;
        MultiGeneSingleOmegaArgParse args(cmd);
        cmd.parse();
        d.chain_driver =
            unique_ptr<D>(new D(cmd.chain_name(),
                                          args.app.every.getValue(),
                                          args.app.until.getValue()));
        d.model = unique_ptr<M>(new M(args.app.alignment.getValue(),
                                      args.app.treefile.getValue(),
                                      myid, nprocs));

        double omegahypermean;
        double omegahyperinvshape;
        d.model->SetAcrossGenesModes(args.blmode(),
                                   args.nucmode(),
                                   args.omegamode(omegahypermean, omegahyperinvshape));
        d.model->SetOmegaHyperParameters(omegahypermean, omegahyperinvshape);
        d.model->Allocate();
        d.model->Update();
        return d;
    }
}

int main(int argc, char *argv[]) {
    int myid = 0;
    int nprocs = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


    ChainCmdLine cmd{argc, argv, "MultiGeneSingleOmega", ' ', "0.1"};

    if(!myid) {
        auto d = load_appdata<ChainDriver, MultiGeneSingleOmegaModelMaster>(cmd, myid, nprocs);
        d.chain_driver->add(*d.model);
        d.chain_driver->go();
    }
    else {
        auto d = load_appdata<SlaveChainDriver, MultiGeneSingleOmegaModelSlave>(cmd, myid, nprocs);
        d.chain_driver->add(*d.model);
        d.chain_driver->go();
    }
    MPI_Finalize();
}
