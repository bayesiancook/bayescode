#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#pragma once

using spdlog::debug;
using spdlog::error;
using spdlog::info;
using spdlog::warn;

#define LIB_INFO(...) SPDLOG_INFO(__VA_ARGS__)
#define LIB_WARN(...) SPDLOG_WARN(__VA_ARGS__)
#define LIB_ERROR(...) SPDLOG_ERROR(__VA_ARGS__)
#define LIB_DEBUG(...) SPDLOG_DEBUG(__VA_ARGS__)
#define LIB_TRACE(...) SPDLOG_TRACE(__VA_ARGS__)
#define LIB_FAIL(...) { SPDLOG_CRITICAL(__VA_ARGS__); exit(1); }