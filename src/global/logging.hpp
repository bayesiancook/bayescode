#pragma once

#ifndef NDEBUG
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#endif

#include <iostream>
#include "spdlog/fmt/ostr.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"
using std::string;

// Type for logger pointers
using logger_t = std::shared_ptr<spdlog::logger>;

// Global variables to simplify formatting
const std::vector<int> color_codes{33, 34, 35, 36, 91, 92, 93, 94, 95, 96};
const string dim_code{"\e[0m\e[2m"}, bold_code{"\e[0m\e[1m"}, normal_code{"\e[0m"};
// const string colored_pattern_prefix{"[" + dim_code + "%Y-%m-%d %H:%M:%S.%e" + normal_code + "] ["
// + dim_code + "%P" + normal_code + "] [" + dim_code + "%n" + normal_code + "] [" + dim_code +
// "%^%l%$" + normal_code + "] "};
const string long_timestamp = "[" + dim_code + "%Y-%m-%d %H:%M:%S.%e" + normal_code + "] ";
const string short_timestamp = "[" + dim_code + "%m-%d %H:%M:%S" + normal_code + "] ";
const string process_block = "[" + dim_code + "%P" + normal_code + "] ";
const string source_block = "[" + dim_code + "%n" + normal_code + "] ";
const string type_block = "[" + dim_code + "%^%l%$" + normal_code + "] ";
const string fileline_block = "[" + dim_code + "%@" + normal_code + "] ";

const string short_prefix = short_timestamp + type_block;

// Factory functions to create loggers
inline logger_t stdout_logger(string name) {
    auto result = spdlog::get(name);
    if (result.get() == nullptr) {
        result = spdlog::stdout_color_mt(name);
        result->set_pattern(short_prefix + source_block + "%v");
    }
    return result;
}

inline logger_t mpi_stdout_logger(string name) {
    auto result = spdlog::get(name);
    if (result.get() == nullptr) {
        result = spdlog::stdout_color_mt(name);
        result->set_pattern(short_prefix + process_block + source_block + "%v");
    }
    return result;
}

inline logger_t file_logger(string name, string filename) {
    auto result = spdlog::get(name);
    if (result.get() == nullptr) {
        result = spdlog::basic_logger_mt(name, filename);
        result->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%P] [%n] [%^%l%$] %v");
    }
    return result;
}

inline logger_t global_logger() {
    auto result = spdlog::get("global");
    if (result == nullptr) {
        result = spdlog::stdout_color_mt("global");
        result->set_pattern(short_prefix + fileline_block + "%v");
    }
    return result;
}

// Macros for global logging
#define INFO(...) SPDLOG_LOGGER_INFO(global_logger(), __VA_ARGS__)
#define WARNING(...) SPDLOG_LOGGER_WARN(global_logger(), __VA_ARGS__)
#define ERROR(...) SPDLOG_LOGGER_ERROR(global_logger(), __VA_ARGS__)
#define DEBUG(...) SPDLOG_LOGGER_DEBUG(global_logger(), __VA_ARGS__)
#define TRACE(...) SPDLOG_LOGGER_TRACE(global_logger(), __VA_ARGS__)
#define FAIL(...)                                             \
    {                                                         \
        SPDLOG_LOGGER_CRITICAL(global_logger(), __VA_ARGS__); \
        exit(1);                                              \
    }
#ifndef NDEBUG
#define assert_warn(condition, ...) \
    if (not(condition)) { WARNING(__VA_ARGS__); }
#else
#define assert_warn(...)
#endif

// utility functions for debug
template <class T>
string vector_to_string(const std::vector<T>& vec) {
    if (vec.size() == 0) {
        return "{}";
    } else {
        string result = "{";
        for (size_t i = 0; i < vec.size() - 1; i++) { result += std::to_string(vec[i]) + ", "; }
        result += std::to_string(vec.back());
        return result + "}";
    }
}
