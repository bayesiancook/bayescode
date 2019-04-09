#include "doctest.h"

#include <sstream>
#include "param_classes.hpp"
#include "parameter_parser.hpp"

TEST_CASE("Basic parser test") {
    std::stringstream ss{"shared"};
    param::Config<double, no_hyper> config;
    ss >> config;
    CHECK(config.mode == param::shared);

    ss.clear();
    ss.str("shrunk");
    config = param::Config<double, no_hyper>();
    ss >> config;
    CHECK(config.mode == param::shrunk);

    ss.clear();
    ss.str("fixed:0.1");
    config = param::Config<double, no_hyper>();
    ss >> config;
    CHECK(config.mode == param::fixed);
    CHECK(config.value == 0.1);

    ss.clear();
    ss.str("independent/0.1");
    param::Config<double, hyper_rate> config2;
    ss >> config2;
    CHECK(config2.mode == param::independent);
    CHECK(config2.hyper.rate == 0.1);
}