#pragma once
#include <string>
#include <vector>
#include "ChainComponent.hpp"
#include "RunToggle.hpp"

class ChainDriver {
    static std::string get_first_token(std::istream& is) {
        std::string result;
        is >> result;
        return result;
    }

  public:
    ChainDriver(std::string name, int every, int until, int size = 0)
        : name(name), toggle(name + ".run"), every(every), until(until), size(size) {}

    void go() {
        if (size == 0)
            for (auto c : components) c->start();
        while (toggle.check() and (until == -1 or size < until)) {
            for (int i = 0; i < every; i++)
                for (auto c : components) c->move(size * every + i);
            for (auto c : components) c->savepoint(size);
            size++;
        }
        for (auto c : components) c->end();
    }

    void add(ChainComponent& component) { components.push_back(&component); }

    void serialize(std::ostream& os) const {
        os << name << "\t" << every << "\t" << until << "\t" << (size + 1);
    }

    ChainDriver(std::istream& is) : name(get_first_token(is)), toggle(name + ".run") {
        is >> every;
        is >> until;
        is >> size;
    }

    static void fake_read(std::istream& is) {
        std::string tmp;
        is >> tmp >> tmp >> tmp >> tmp;
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
