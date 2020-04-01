#include "MultiGeneChain.hpp"
#include <fstream>
#include <iostream>
#include "Chrono.hpp"
#include "MultiGeneProbModel.hpp"
using namespace std;

// c++11
#define nullptr 0

MultiGeneChain::MultiGeneChain(int inmyid, int innprocs)
    : Chain(), myid(inmyid), nprocs(innprocs), maxtime(0) {}

void MultiGeneChain::SavePoint() {
    if (saveall) {
        if (!myid) {
            ofstream chain_os((name + ".chain").c_str(), ios_base::app);
            GetMultiGeneModel()->MasterToStream(chain_os);
        } else {
            GetMultiGeneModel()->SlaveToStream();
        }
    }
    size++;
}

void MultiGeneChain::Reset(int force) {
    size = 0;
    if (!myid) {
        MakeFiles(force);
    }
    Save();
}

void MultiGeneChain::MakeFiles(int force) {
    Chain::MakeFiles(force);
    ofstream nameos((name + ".genelist").c_str());
    GetMultiGeneModel()->PrintGeneList(nameos);
    nameos.close();
}

void MultiGeneChain::Move() {
    for (int i = 0; i < every; i++) {
        GetMultiGeneModel()->Move();
    }
    SavePoint();
    Save();
    if (!myid) {
        Monitor();
    }
}

void MultiGeneChain::Start() {
    if (!myid) {
        ofstream run_os((name + ".run").c_str());
        run_os << 1 << '\n';
        run_os.close();
        if (maxtime)    {
            global_chrono.Start();
        }
    }
    Run();
}

void MultiGeneChain::MasterSendRunningStatus(int status) {
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

int MultiGeneChain::SlaveReceiveRunningStatus() {
    int status;
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return status;
}

int MultiGeneChain::GetRunningStatus()  {
    int ret = Chain::GetRunningStatus();
    if (maxtime)    {
        global_chrono.Stop();
        ret &= ((global_chrono.GetTime() / 3600000) < maxtime);
        global_chrono.Start();
    }
    return ret;
}

void MultiGeneChain::Run() {
    if (!myid) {
        while ((GetRunningStatus() != 0) && ((until == -1) || (size <= until))) {
            MasterSendRunningStatus(1);
            Chrono chrono;
            chrono.Start();
            Move();
            chrono.Stop();

            ofstream check_os((name + ".time").c_str());
            check_os << chrono.GetTime() << '\t' << global_chrono.GetTime() << '\n';
            cerr << global_chrono.GetTime() / 3600000 << '\t' << maxtime << '\n';
        }
        MasterSendRunningStatus(0);
        ofstream run_os((name + ".run").c_str());
        run_os << 0 << '\n';
    } else {
        while (SlaveReceiveRunningStatus()) {
            Move();
        }
    }
}
