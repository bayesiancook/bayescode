#pragma once

#include <functional>
#include <string>
#include "ChainComponent.hpp"
#include "ChainDriver.hpp"

class ChainCheckpoint : public ChainComponent {
    std::string filename;
    std::function<void(std::ostream &)> serialize_model;
    ChainDriver &cd;

  public:
    template <class T>
    ChainCheckpoint(std::string filename, ChainDriver &cd, T &model)
        : filename(filename),
          serialize_model([&model](std::ostream &os) { model.ToStream(os); }),
          cd(cd) {}

    void savepoint(int) override {
        std::ofstream os{filename};
        cd.serialize(os);
        os << "\n";
        serialize_model(os);
    }
};
