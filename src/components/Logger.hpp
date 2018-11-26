#pragma once

#include <functional>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>
#include "utils/date.h"

class Token {
    std::function<std::string()> raw_string;
    std::function<std::string(std::string)> modified_string{[](std::string s) { return s; }};

  public:
    template <class... Args>
    Token(std::set<int> modifiers, std::string format, Args&&... args) {
        // formatting to intermediate string
        char* buf = nullptr;
        int res = asprintf(&buf, format.c_str(), std::forward<Args>(args)...);
        if (res == -1) {
            fprintf(stderr, "Error in Logger::message: something went wrong in asprintf!\n");
        }

        // storing raw string in function
        std::string buf_string(buf);
        free(buf);
        raw_string = [buf_string]() { return buf_string; };

        std::string modif_prefix, modif_suffix;
        if (modifiers.size() > 0) {
            for (auto modifier : modifiers) {
                modif_prefix += "\e[" + std::to_string(modifier) + "m";
            }
            modif_suffix += "\e[0m";
        }
        modified_string = [modif_prefix, modif_suffix](
            std::string s) { return modif_prefix + s + modif_suffix; };
    }

    template <class F>
    Token(F f) : raw_string(f) {}
    Token(std::string s) : raw_string([s]() { return s; }) {}
    Token(const char* s) : raw_string([s]() { return std::string(s); }) {}

    std::string str() const { return modified_string(raw_string()); }

    size_t size() const { return raw_string().size(); }
};

std::string timestamp_time() {
    return date::format("%H:%M:%S",
        std::chrono::time_point_cast<std::chrono::seconds>(std::chrono::system_clock::now()));
}

std::string timestamp_date() { return date::format("%Y %b %d", std::chrono::system_clock::now()); }

std::string timestamp_full() {
    return date::format("%Y-%m-%dT%T%z", std::chrono::system_clock::now());
}

// quick functions for specific token types
template <class... Args>
Token bold(std::string format, Args&&... args) {
    return Token({1}, format, std::forward<Args>(args)...);
}

template <class... Args>
Token red(std::string format, Args&&... args) {
    return Token({31}, format, std::forward<Args>(args)...);
}

template <class... Args>
Token bold_red(std::string format, Args&&... args) {
    return Token({1, 31}, format, std::forward<Args>(args)...);
}

class MessageFormat {
    std::vector<Token> _prefix, _suffix, _line_prefix;

  public:
    MessageFormat(std::vector<Token> prefix, std::vector<Token> suffix = {},
        std::vector<Token> line_prefix = {})
        : _prefix(prefix), _suffix(suffix), _line_prefix(line_prefix) {}

    std::string prefix() const {
        return std::accumulate(_prefix.begin(), _prefix.end(), std::string(),
            [](std::string acc, Token t) { return acc + t.str(); });
    }
    std::string suffix() const {
        return std::accumulate(_suffix.begin(), _suffix.end(), std::string(),
            [](std::string acc, Token t) { return acc + t.str(); });
    }
    std::string line_prefix() const {
        return std::accumulate(_line_prefix.begin(), _line_prefix.end(), std::string(),
            [](std::string acc, Token t) { return acc + t.str(); });
    }
};

class Logger {
    std::ostream& os;

  public:
    Logger(std::ostream& os) : os(os) {}

    template <class... Args>
    void print(const std::string& format, Args&&... args) const {
        std::cout << Token({}, format, std::forward<Args>(args)...).str();
    }

    template <class... Args>
    void message(MessageFormat& message_format, std::string format, Args&&... args) const {
        // replace all \n with \n + line_prefix
        auto message = Token({}, format, std::forward<Args>(args)...).str();
        auto line_prefix = "\n" + message_format.line_prefix();
        size_t pos = 0;
        while (true) {
            pos = message.find('\n', pos + 1);
            if (pos == std::string::npos) { break; }
            message.replace(pos, 1, line_prefix);
        }

        std::cout << message_format.prefix() << message << message_format.suffix();
    }
};