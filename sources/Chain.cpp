#include "Chain.hpp"
#include <fstream>
#include <iostream>
#include "Chrono.hpp"
#include "ProbModel.hpp"
using namespace std;

// c++11
#define nullptr 0

Chain::Chain() {
    every = 1;
    until = -1;
    saveall = 1;
    size = 0;
    model = nullptr;
    name = "";
}

void Chain::MakeFiles(int force) {
    if (ifstream((name + ".param").c_str()) && (force == 0)) {
        cerr << "already existing chain, cannot override (unless in forcing mode)\n";
        exit(1);
    }
    ofstream param_os((name + ".param").c_str());
    if (saveall) { ofstream chain_os((name + ".chain").c_str()); }
    ofstream mon_os((name + ".monitor").c_str());
    ofstream trace_os((name + ".trace").c_str());
    model->TraceHeader(trace_os);
}

void Chain::Monitor() {
    ofstream trace_os((name + ".trace").c_str(), ios_base::app);
    model->Trace(trace_os);
    ofstream mon_os((name + ".monitor").c_str());
    ofstream mon_det_os((name + ".details").c_str());
    model->Monitor(mon_os);
}

void Chain::SavePoint() {
    if (saveall) {
        ofstream chain_os((name + ".chain").c_str(), ios_base::app);
        model->ToStream(chain_os);
    }
    size++;
}

void Chain::Reset(int force) {
    size = 0;
    MakeFiles(force);
    Save();
}

void Chain::Move() {
    for (int i = 0; i < every; i++) { model->Move(); }
    SavePoint();
    Save();
    Monitor();
}

void Chain::Start() {
    ofstream run_os((name + ".run").c_str());
    run_os << 1 << '\n';
    run_os.close();
    Run();
}

int Chain::GetRunningStatus() {
    ifstream run_is((name + ".run").c_str());
    int run;
    run_is >> run;
    return run;
}

void Chain::Run() {
    while ((GetRunningStatus() != 0) && ((until == -1) || (size <= until))) {
        Chrono chrono;
        chrono.Start();
        Move();
        chrono.Stop();
        ofstream check_os((name + ".time").c_str());
        check_os << chrono.GetTime() << '\n';
    }
    ofstream run_os((name + ".run").c_str());
    run_os << 0 << '\n';
}
