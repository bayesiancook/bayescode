#include "Chain.hpp"
#include <fstream>
#include "Chrono.hpp"
#include "ProbModel.hpp"
using namespace std;

void ChainDriver::go() {
    chainfile->write_header();
    tracefile->write_header();

    lifecycle_handler->Init();

    tracefile->write_line();
    chainfile->write_line();

    Run();

    lifecycle_handler->End();
}

void ChainDriver::Move() {
    for (int i = 0; i < every; i++) {
        model->Move();
    }
    tracefile->write_line();
    chainfile->write_line();
    monitorfile->write_line();

    size++;

    lifecycle_handler->EndMove();
}

bool ChainDriver::IsRunning() { return run_toggle->check(); }

void ChainDriver::Run() {
    while (IsRunning() && ((until == -1) || (size <= until))) {
        Chrono chrono;
        chrono.Reset();
        chrono.Start();
        Move();
        chrono.Stop();
    }
}
