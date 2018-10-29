#pragma once

#include <mpi.h>
#include <memory>
#include <vector>

/*==================================================================================================
  Process
  Class representing a MPI process. Provides size, rank and messaging functionality.
==================================================================================================*/
class Process {
    static const std::vector<int> colors;
    static const std::string bold, normal;
    const std::string color, prefix;

  public:
    Process(int rank, int size)
        : color("\e[0m\e[" + std::to_string(colors[rank % colors.size()]) + "m"),
          prefix{bold + "[" + color + "%d" + bold + "/" + color + "%d" + bold + "] " + normal},
          rank(rank),
          size(size) {}

    const int rank, size;

    template <class... Args>
    void message(const std::string& format, Args&&... args) {
        std::string format2 = prefix + format + "\n";
        printf(format2.c_str(), rank, size, std::forward<Args>(args)...);
    }
};

const std::vector<int> Process::colors{31, 32, 33, 34, 35, 36, 91, 92, 93, 94, 95, 96};
const std::string Process::bold{"\e[0m\e[1m"}, Process::normal{"\e[0m"};

/*==================================================================================================
  MPI::p
  Global Process object.
==================================================================================================*/
namespace MPI {
    std::unique_ptr<Process> p{nullptr};
};

/*==================================================================================================
  mpi_run
  Wrapper around a main-like function that initializes MPI and sets MPI::p to correspond to the
  local MPI process
==================================================================================================*/
template <class F, class... Args>
void mpi_run(int argc, char** argv, F f) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI::p = std::make_unique<Process>(rank, size);
    MPI::p->message("Started MPI process");
    f(argc, argv);
    MPI::p->message("End of MPI process");
    MPI_Finalize();
}