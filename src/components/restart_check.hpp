#pragma once

#include "components/Tracer.hpp"

template <class Model>
bool check_restart(Model& model, std::string tracefile) {
    // get stats from new model
    Tracer tracer(model, &Model::declare_stats);
    std::stringstream ss;
    tracer.write_line(ss);

    // get stats from trace file
    std::ifstream trace{tracefile};
    std::string line, nonempty_line;
    while (std::getline(trace, line)) {
        if (line != "\n") { nonempty_line = line; }
    }

    if (("\n" + nonempty_line) != ss.str()) {
        std::cerr << "Warning: statistics from restarted model do not match statistics from "
                     "tracefile "
                  << tracefile << "\n";
        std::cerr << "Restarted statistics are:\n";
        tracer.write_header(std::cerr);
        tracer.write_line(std::cerr);
        std::cerr << "\nStatistics loaded from trace:\n";
        tracer.write_header(std::cerr);
        std::cerr << "\n" << nonempty_line << "\n";
        return false;
    } else {
        std::cout << "Restart check passed successfully!\n";
        std::cout << "Statistics are:\n";
        tracer.write_header(std::cout);
        tracer.write_line(std::cout);
        std::cerr << "\n";
        return true;
    }
}