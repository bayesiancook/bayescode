#include "doctest.h"

#include "parameter_parser.hpp"

TEST_CASE("Basic parser test") {
    std::stringstream ss{"test:0.1/0.2;0.3:yolo"};
    MultiGeneParameter<double, no_hyper> test;
    ss >> test;
    CHECK(false);
}