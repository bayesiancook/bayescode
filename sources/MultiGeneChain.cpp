#include "MultiGeneChain.hpp"
#include <ctime>
#include <fstream>
#include <iostream>
#include <list>
#include "Chrono.hpp"
#include "MultiGeneProbModel.hpp"
#include "monitoring.hpp"

std::unique_ptr<MonitorManager> gm(new MonitorManager());

using namespace std;

// c++11
#define nullptr 0

MultiGeneChain::MultiGeneChain(int inmyid, int innprocs)
    : Chain(), myid(inmyid), nprocs(innprocs) {}

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

void MultiGeneChain::Run() {
    ofstream mvfile(name + "_p" + to_string(myid) + "_" + to_string(size) + "to" +
                    to_string(until) + ".moveaccept");
    ofstream mvfiletime(name + "_p" + to_string(myid) + "_" + to_string(size) + "to" +
                        to_string(until) + ".movetime");
    ofstream mvfiletottime(name + "_p" + to_string(myid) + "_" + to_string(size) + "to" +
                           to_string(until) + ".movetottime");

    int first_iteration = size + 1;

    auto write_line = [this, &mvfile, &mvfiletime, &mvfiletottime]() {
        stringstream line, timeline, tottimeline;
        for (auto& monitor : gm->monitors) {
            if (line.str() != "") {
                line << '\t';
                timeline << '\t';
                tottimeline << '\t';
            }
            line << dynamic_cast<MeanMonitor<double>*>(monitor.second.get())->tmp_mean();
            timeline << dynamic_cast<MeanMonitor<double>*>(monitor.second.get())->tmp_time_mean();
            tottimeline
                << dynamic_cast<MeanMonitor<double>*>(monitor.second.get())->tmp_time_total();
            dynamic_cast<MeanMonitor<double>*>(monitor.second.get())->tmp_reset();
        }
        mvfile << line.str() << "\n";
        mvfiletime << timeline.str() << "\n";
        mvfiletottime << tottimeline.str() << "\n";
    };

    auto write_header = [this, &mvfile, &mvfiletime, &mvfiletottime]() {
        stringstream header;
        for (auto& monitor : gm->monitors) {
            if (header.str() != "") {
                header << '\t';
            }
            header << monitor.first;
        }
        mvfile << header.str() << "\n";
        mvfiletime << header.str() << "\n";
        mvfiletottime << header.str() << "\n";
    };

    list<double> last_it_times;
    time_t max_time = 24 * 60 * 60;  // 24 hours * 60 mintues * 60 seconds
    time_t start_time = std::time(nullptr);

    auto mean_and_trim = [this, &last_it_times](double it_time) -> double {
        last_it_times.push_back(it_time);
        double sum = 0;
        int count = 0;
        for (auto t : last_it_times) {
            count++;
            sum += t;
        }
        if (last_it_times.size() > 5) {
            last_it_times.pop_front();
        }
        return sum / count;
    };

    //=============================================================================================
    //  MAIN LOOP
    //=============================================================================================
    if (!myid) {
        std::cout << "Computation started at " << std::asctime(std::localtime(&start_time));
        while ((GetRunningStatus() != 0) && ((until == -1) || (size <= until))) {
            MasterSendRunningStatus(1);
            Chrono chrono;
            chrono.Start();
            Move();
            chrono.Stop();

            double it_time = chrono.GetTime();
            double mean_it_time = mean_and_trim(it_time) / 1000;
            time_t current_time = std::time(nullptr);
            time_t total_time = current_time - start_time;
            int remaining_its = floor((max_time - total_time) / mean_it_time);

            ofstream check_os((name + ".time").c_str());
            std::string time_str = std::asctime(std::localtime(&current_time));
            time_str.pop_back();
            check_os << it_time << '\n';
            cerr << "[" << time_str << "] Iteration " << size - 1 << ": " << it_time / 1000
                 << "s (mean it time: " << mean_it_time << "s; predicting " << remaining_its - 2
                 << " more iterations)\n";

            if (size == first_iteration) {
                write_header();
            }
            write_line();

            if (total_time + 3 * mean_it_time > max_time) {
                printf("Appraoching max time! Stopping computation.\n");
                break;
            }
        }
        MasterSendRunningStatus(0);
        ofstream run_os((name + ".run").c_str());
        run_os << 0 << '\n';
    } else {
        while (SlaveReceiveRunningStatus()) {
            Move();
            write_line();
        }
    }
}
