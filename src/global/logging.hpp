#pragma once

#include <iostream>
#include "spdlog/fmt/ostr.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

using std::string;

using logger_t = std::shared_ptr<spdlog::logger>;

logger_t stdout_logger(string name);

logger_t file_logger(string name, string filename);

extern logger_t global_logger;
extern const std::vector<int> color_codes;
extern const std::string dim_code, bold_code, normal_code;

#define INFO(...) SPDLOG_LOGGER_INFO(global_logger, __VA_ARGS__)
#define WARNING(...) SPDLOG_LOGGER_WARN(global_logger, __VA_ARGS__)
#define ERROR(...) SPDLOG_LOGGER_ERROR(global_logger, __VA_ARGS__)
#define DEBUG(...) SPDLOG_LOGGER_DEBUG(global_logger, __VA_ARGS__)
#define TRACE(...) SPDLOG_LOGGER_TRACE(global_logger, __VA_ARGS__)
#define FAIL(...)                     \
    {                                 \
        SPDLOG_CRITICAL(__VA_ARGS__); \
        exit(1);                      \
    }