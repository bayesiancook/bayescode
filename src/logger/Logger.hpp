#pragma once

#include <iostream>
#include <numeric>
#include <set>
#include <vector>

class Token {
    size_t _size;
    std::string _value;

  public:
    template <class... Args>
    Token(std::set<int> modifiers, std::string format, Args&&... args) {
        // formatting to intermediate string
        char* buf = nullptr;
        int res = asprintf(&buf, format.c_str(), std::forward<Args>(args)...);
        if (res == -1) {
            fprintf(stderr, "Error in Logger::message: something went wrong in asprintf!\n");
        }
        _value = buf;
        free(buf);
        _size = _value.size();  // size BEFORE modifiers (ie, in terms of displayed characters)

        if (modifiers.size() > 0) {
            for (auto modifier : modifiers) {
                _value = "\e[" + std::to_string(modifier) + "m" + _value;
            }
            _value += "\e[0m";
        }
    }

    std::string str() const { return _value; }

    size_t size() const { return _size; }
};

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
  public:
    template <class... Args>
    void print(const std::string& format, Args&&... args) const {
        std::cout << Token({}, format, std::forward<Args>(args)...).str();
    }

    template <class... Args>
    void message(MessageFormat& message_format, std::string format, Args&&... args) const {
        std::cout << message_format.prefix() << Token({}, format, std::forward<Args>(args)...).str()
                  << message_format.suffix();
    }
};