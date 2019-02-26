#pragma once

#include <iostream>
#include "spdlog/fmt/ostr.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

using std::string;

using logger_t = std::shared_ptr<spdlog::logger>;

inline logger_t stdout_logger(string name) { return spdlog::stdout_color_mt(name); }

inline logger_t file_logger(string name, string filename) {
    return spdlog::basic_logger_mt(name, filename);
}

extern logger_t lib_logger;

#define INFO(...) SPDLOG_LOGGER_INFO(lib_logger, __VA_ARGS__)
#define WARN(...) SPDLOG_LOGGER_WARN(lib_logger, __VA_ARGS__)
#define ERROR(...) SPDLOG_LOGGER_ERROR(lib_logger, __VA_ARGS__)
#define DEBUG(...) SPDLOG_LOGGER_DEBUG(lib_logger, __VA_ARGS__)
#define TRACE(...) SPDLOG_LOGGER_TRACE(lib_logger, __VA_ARGS__)
#define FAIL(...)                     \
    {                                 \
        SPDLOG_CRITICAL(__VA_ARGS__); \
        exit(1);                      \
    }