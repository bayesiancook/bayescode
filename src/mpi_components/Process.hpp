#pragma once

#include <mpi.h>
#include <vector>

/*==================================================================================================
  Process
  Class representing a MPI process. Provides size, rank and messaging functionality.
==================================================================================================*/
class Process {
    static const std::vector<int> colors;

  public:
    int rank{-1}, size{-1};

    template <class... Args>
    void message(const std::string& format, Args&&... args) {
        std::string color = "\e[0m\e[" + std::to_string(colors[rank % colors.size()]) + "m";
        std::string bold = "\e[0m\e[1m";
        std::string normal = "\e[0m";
        std::string format2 = bold + "[" + color + "%d" + bold + "/" + color + "%d" + bold + "] " +
                              normal + format + "\n";
        printf(format2.c_str(), rank, size, std::forward<Args>(args)...);
    }
};

const std::vector<int> Process::colors{31, 32, 33, 34, 35, 36, 91, 92, 93, 94, 95, 96};

namespace MPI {
    Process p;
};

template <class F, class... Args>
void mpi_run(int argc, char** argv, F f) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI::p.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI::p.size);
    MPI::p.message("Started MPI process");
    f(argc, argv);
    MPI::p.message("End of MPI process");
    MPI_Finalize();
}