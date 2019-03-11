#pragma once

#include "ChainComponent.hpp"
#include "global/logging.hpp"

class ConsoleLogger : public ChainComponent {
    logger_t logger{stdout_logger("chain_logger")};

  public:
    void start() override { logger->info("Started"); }
    void move(int i) override { logger->info("Move {}", i); }
    void savepoint(int i) override { logger->info("Savepoint {}", i); }
    void end() override { logger->info("Ended"); }
};
