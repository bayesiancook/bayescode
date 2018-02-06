#include "MultiGeneChain.hpp"
#include <fstream>
#include <iostream>
#include "MultiGeneProbModel.hpp"
#include "Chrono.hpp"
using namespace std;

// c++11
#define nullptr 0 

MultiGeneChain::MultiGeneChain(int inmyid, int innprocs) : Chain(), myid(inmyid), nprocs(innprocs) {}

void MultiGeneChain::SavePoint() {
    if (! myid) {
        ofstream chain_os((name + ".chain").c_str(), ios_base::app);
        GetMultiGeneModel()->MasterToStream(chain_os);
    }
    else    {
        GetMultiGeneModel()->SlaveToStream();
    }
    size++;
}

void MultiGeneChain::Reset(int force)   {
    size = 0;
    if (! myid) {
        MakeFiles(force);
    }
    Save();
}

void MultiGeneChain::Move() {

    for (int i = 0; i < every; i++) {
        GetMultiGeneModel()->Move();
    }
    SavePoint();
    Save();
    if (!myid)  {
        Monitor();
    }
}

void MultiGeneChain::Start() {
    if (! myid) {
        ofstream run_os((name + ".run").c_str());
        run_os << 1 << '\n';
        run_os.close();
    }
    Run();
}

void MultiGeneChain::MasterSendRunningStatus(int status)  {
    MPI_Bcast(&status,1,MPI_INT,0,MPI_COMM_WORLD);
}

int MultiGeneChain::SlaveReceiveRunningStatus() {
    int status;
    MPI_Bcast(&status,1,MPI_INT,0,MPI_COMM_WORLD);
    return status;
}

void MultiGeneChain::Run() {

    if (! myid) {

        int i = 0;
        MeasureTime timer;

        while ((GetRunningStatus() != 0) && ((until == -1) || (size <= until))) {

            MasterSendRunningStatus(1);
            Chrono chrono;
            chrono.Reset();
            chrono.Start();
            Move();
            chrono.Stop();

            timer << "Iteration " << i * every << ". ";
            timer.print<0>();
            i++;

            /*
            ofstream check_os((name + ".time").c_str());
            check_os << chrono.GetTime() / 1000 << '\n';
            */
        }
        MasterSendRunningStatus(0);
        ofstream run_os((name + ".run").c_str());
        run_os << 0 << '\n';
    }
    else    {

        while (SlaveReceiveRunningStatus()) {
            Move();
        }
    }
}
