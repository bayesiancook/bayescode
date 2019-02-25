#pragma once

#include "ChainComponent.hpp"
#include "Logger.hpp"

class ConsoleLogger : public ChainComponent {
  public:
    void start() override { info("Started"); }
    void move(int i) override { info("Move {}", i); }
    void savepoint(int i) override { info("Savepoint {}", i); }
    void end() override { info("Ended"); }
};
