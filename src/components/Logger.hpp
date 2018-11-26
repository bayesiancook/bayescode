#pragma once

#include <iostream>
#include <numeric>
#include <set>
#include <vector>

class Token {
    std::string _value;
    size_t _size;

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

    Token(std::string s) : _value(s), _size(_value.size()) {}
    Token(const char* s) : _value(s), _size(_value.size()) {}

    std::string str() const { return _value; }

    size_t size() const { return _size; }
};

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