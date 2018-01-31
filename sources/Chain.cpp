#include "Chain.hpp"
#include <fstream>
#include <iostream>
#include "Chrono.hpp"
#include "ProbModel.hpp"
using namespace std;

void Chain::MakeFiles(bool force) {
    if (ifstream((name + ".param").c_str()) && !force) {
        cerr << "already existing chain, cannot override (unless in forcing mode)\n";
        exit(1);
    }
    ofstream param_os((name + ".param").c_str());
    ofstream chain_os((name + ".chain").c_str());
    ofstream mon_os((name + ".monitor").c_str());
    ofstream trace_os((name + ".trace").c_str());
    model->TraceHeader(trace_os);
}

void Chain::Monitor() {
    ofstream trace_os((name + ".trace").c_str(), ios_base::app);
    model->Trace(trace_os);
    ofstream mon_os((name + ".monitor").c_str());
    model->Monitor(mon_os);
}

void Chain::SavePoint() {
    ofstream chain_os((name + ".chain").c_str(), ios_base::app);
    model->ToStream(chain_os);
    size++;
}

void Chain::Reset(bool force) {
    size = 0;
    MakeFiles(force);
    // Save();
}

void Chain::Move() {
    for (int i = 0; i < every; i++) {
        model->Move();
    }
    SavePoint();
    // Save();
    Monitor();
}

bool Chain::IsRunning() { return run_toggle->check(); }

void Chain::Run() {
    while (IsRunning() && ((until == -1) || (size <= until))) {
        Chrono chrono;
        chrono.Reset();
        chrono.Start();
        Move();
        chrono.Stop();
    }
}
