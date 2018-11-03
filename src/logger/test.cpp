#include "Logger.hpp"

int main() {
    Logger l;
    l.print(std::cout, "Hello %d!\n", 13);
}