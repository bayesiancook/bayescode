#pragma once

#include <iostream>
#include <set>

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


class Logger {
  public:
    template <class... Args>
    void print(std::ostream& os, const std::string& format, Args&&... args) {
        os << Token({1, 31}, format, std::forward<Args>(args)...).str();
    }
};