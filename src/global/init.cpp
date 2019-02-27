#include "init.hpp"

// Global decls
logger_t global_logger = stdout_logger("global");
const std::vector<int> color_codes{33, 34, 35, 36, 91, 92, 93, 94, 95, 96};
const string dim_code{"\e[0m\e[2m"}, bold_code{"\e[0m\e[1m"}, normal_code{"\e[0m"};
const string colored_pattern_prefix =
    "[" + dim_code + "%Y-%m-%d %H:%M:%S.%e" + normal_code + "] [" + dim_code + "%P" + normal_code +
    "] [" + dim_code + "%n" + normal_code + "] [" + dim_code + "%^%l%$" + normal_code + "] ";

// Formatting-related functions to keep it all in one place

logger_t stdout_logger(string name) {
    auto result = spdlog::get(name);
    if (result.get() == nullptr) {
        result = spdlog::stdout_color_mt(name);
        result->set_pattern(colored_pattern_prefix + "%v");
    }
    return result;
}

logger_t file_logger(string name, string filename) {
    auto result = spdlog::get(name);
    if (result.get() == nullptr) {
        result = spdlog::basic_logger_mt(name, filename);
        result->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%P] [%n] [%^%l%$] %v");
    }
    return result;
}

// All things to be run before the main
// (typically to initialize global library things)
class init_before_main {
  public:
    init_before_main() {
        global_logger->set_pattern(colored_pattern_prefix + "[" + dim_code + "%@" + normal_code +
                                   "] " + normal_code + "%v");

        // random seed
        Random::InitRandom();
        INFO("Seed: {}", Random::GetSeed());
    }
};

static init_before_main init;