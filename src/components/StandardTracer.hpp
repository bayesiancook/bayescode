#pragma once
#include <fstream>
#include "ChainComponent.hpp"
#include "Tracer.hpp"

class StandardTracer : public ChainComponent {
    Tracer model_tracer;
    Tracer stats_tracer;
    std::string chain_name;

  public:
    template <class M>
    StandardTracer(M& m, std::string chain_name)
        : model_tracer(m), stats_tracer(m, processing::HasTag<Stat>()), chain_name(chain_name) {}

    void start() final {
        std::ofstream model_os{chain_file(chain_name), std::ios_base::trunc};
        model_tracer.write_header(model_os);
        std::ofstream stats_os{chain_name + ".trace", std::ios_base::trunc};
        stats_tracer.write_header(stats_os);
    }

    void savepoint(int) final {
        std::ofstream model_os{chain_file(chain_name), std::ios_base::app};
        model_tracer.write_line(model_os);
        std::ofstream trace_os{chain_name + ".trace", std::ios_base::app};
        stats_tracer.write_line(trace_os);
    }

    static std::string chain_file(std::string chain_name) { return chain_name + ".chain"; }
};
