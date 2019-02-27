#pragma once

#include <memory>
#include <vector>
#include "global/logging.hpp"

/*==================================================================================================
  Process
  Class representing a MPI process. Provides size, rank and messaging functionality.
==================================================================================================*/
class Process {
    const std::string color;
    logger_t logger;

  public:
    Process(int rank, int size)
        : color("\e[0m\e[" + std::to_string(color_codes[rank % color_codes.size()]) + "m"),
          logger(stdout_logger(color + "process-" + std::to_string(rank) + "/" +
                               std::to_string(size) + normal_code)),
          rank(rank),
          size(size) {}

    const int rank, size;

    template <class... Args>
    void message(Args&&... args) {
        logger->info(std::forward<Args>(args)...);
    }
};

/*==================================================================================================
  MPI::p
  Global Process object.
==================================================================================================*/
namespace MPI {
    static std::unique_ptr<Process> p;
};
