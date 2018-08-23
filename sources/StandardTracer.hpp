#include <fstream>
#include "ChainComponent.hpp"

class StandardTracer : public ChainComponent {
    Tracer model_tracer;
    Tracer stats_tracer;
    std::string chain_name;

public:
    template<class M>
    StandardTracer(M& m, std::string chain_name) : model_tracer(m,&M::declare_model),
                                                   stats_tracer(m,&M::declare_stats),
                                                   chain_name(chain_name) {}
    void start() {
        std::ofstream model_os{chain_name + ".chain", std::ios_base::trunc};
        model_tracer.write_header(model_os);
        std::ofstream stats_os{chain_name + ".trace", std::ios_base::trunc};
        stats_tracer.write_header(stats_os);
    }
    void savepoint(int) {
        std::ofstream model_os{chain_name + ".chain", std::ios_base::app};
        model_tracer.write_line(model_os);
        std::ofstream trace_os{chain_name + ".trace", std::ios_base::app};
        stats_tracer.write_line(trace_os);
    }
};
