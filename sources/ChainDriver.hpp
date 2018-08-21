#pragma once
#include <string>
#include <vector>
#include "ChainComponent.hpp"
#include "RunToggle.hpp"

class ChainDriver {
  public:
    ChainDriver(std::string name) : name(name), toggle(name + ".run") {}
    void go() {
        for (auto c : components) c->start();
        while (toggle.check() and (until == -1 or size < until)) {
            for (auto c : components) c->move(size);
            for (auto c : components) c->after_move(size);
            size++;
        }
        for (auto c : components) c->end();
    }

  private:
    std::string name;
    std::vector<ChainComponent*> components;
    RunToggle toggle;
    //! saving frequency (i.e. number of move cycles performed between each point
    //! saved to file)
    int every{1};
    //! intended final size of the chain (until==-1 means no a priori specified
    //! upper limit)
    int until{-1};
    //! current size (number of points saved to file)
    int size{0};
};
