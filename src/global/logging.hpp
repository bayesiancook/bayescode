#pragma once

#include <iostream>
#include "spdlog/spdlog.h"
#include "spdlog/fmt/fmt.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

using std::string;

// Type for logger pointers
using logger_t = std::shared_ptr<spdlog::logger>;

// Global variables to simplify formatting
const std::vector<int> color_codes{33, 34, 35, 36, 91, 92, 93, 94, 95, 96};
const string dim_code{"\e[0m\e[2m"}, bold_code{"\e[0m\e[1m"}, normal_code{"\e[0m"};
const string colored_pattern_prefix{"[" + dim_code + "%Y-%m-%d %H:%M:%S.%e" + normal_code + "] [" +
                                    dim_code + "%P" + normal_code + "] [" + dim_code + "%n" +
                                    normal_code + "] [" + dim_code + "%^%l%$" + normal_code + "] "};

// Factory functions to create loggers
inline logger_t stdout_logger(const string& name) {
    auto result = spdlog::get(name);
    if (result == nullptr) {
        result = spdlog::stdout_color_mt(name);
        result->set_pattern(colored_pattern_prefix + "%v");
    }
    return result;
}

inline logger_t file_logger(const string& name, const string& filename) {
    auto result = spdlog::get(name);
    if (result == nullptr) {
        result = spdlog::basic_logger_mt(name, filename);
        result->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%P] [%n] [%^%l%$] %v");
    }
    return result;
}

inline logger_t global_logger() {
    auto result = spdlog::get("global");
    if (result == nullptr) {
        result = spdlog::stdout_color_mt("global");
        result->set_pattern(colored_pattern_prefix + "[" + dim_code + "%@" + normal_code + "] " +
                            normal_code + "%v");
    }
    return result;
}

// Macros for global logging
#define INFO(...) SPDLOG_LOGGER_INFO(global_logger(), __VA_ARGS__)
#define WARNING(...) SPDLOG_LOGGER_WARN(global_logger(), __VA_ARGS__)
#define ERROR(...) SPDLOG_LOGGER_ERROR(global_logger(), __VA_ARGS__)
#define DEBUG(...) SPDLOG_LOGGER_DEBUG(global_logger(), __VA_ARGS__)
#define TRACE(...) SPDLOG_LOGGER_TRACE(global_logger(), __VA_ARGS__)
#define FAIL(...)                     \
    {                                 \
        SPDLOG_CRITICAL(__VA_ARGS__); \
        exit(1);                      \
    }
