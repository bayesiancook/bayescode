#pragma once

#include <unordered_map>
#include "components/ChainComponent.hpp"
#include "global/logging.hpp"

class OnlineMean {
    double sum{0};
    size_t n{0};

  public:
    void add(double x) {
        sum += x;
        n++;
    }
    double mean() const { return sum / double(n); }
};

class MoveStatsRegistry : public ChainComponent {
    std::unordered_map<std::string, OnlineMean> map;

  public:
    MoveStatsRegistry() = default;
    MoveStatsRegistry(const MoveStatsRegistry&) = delete;

    void add(std::string name, double x) { map[name].add(x); }

    void savepoint(int) final {
        for (auto e : map) {
            DEBUG("Move success: {:>15} -> {:.2f}\%", e.first, 100. * e.second.mean());
        }
    }
};

class MoveReporter {
    std::string name{""};
    MoveStatsRegistry& registry;

  public:
    MoveReporter(std::string name, MoveStatsRegistry& registry) : name(name), registry(registry) {}

    void report(double x) { registry.add(name, x); }
};

class NoReport {
  public:
    void report(double) {}
};