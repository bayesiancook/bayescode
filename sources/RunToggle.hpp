#pragma once
#include <fstream>
#include <string>

class RunToggle {
    std::string filename;

  public:
    RunToggle(std::string filename) : filename(filename) {
        std::ofstream fs;
        fs.open(filename, std::ios_base::in | std::ios_base::out | std::ios_base::trunc);
        fs << 1;
    }

    bool check() {
        std::ifstream fs;
        fs.open(filename, std::ios_base::in);
        return static_cast<char>(fs.peek()) == '1';
    }
};
