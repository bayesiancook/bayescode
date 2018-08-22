#pragma once

#include <iostream>
#include "ChainComponent.hpp"

class ConsoleLogger : public ChainComponent {
  public:
    void start() override { std::cout << "Started\n"; }
    void move(int i) override { std::cout << "Move " << i << "\n"; }
    void savepoint(int i) override { std::cout << "Savepoint " << i << "\n"; }
    void end() override { std::cout << "Ended\n"; }
};
