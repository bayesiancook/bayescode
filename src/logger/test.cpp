#include "Logger.hpp"

int main() {
    MessageFormat error(
        {Token({1}, "["), Token({1, 31}, "error"), Token({1}, "]: ")}, {Token({}, "\n")});

    Logger l;
    l.print("Hello %d!\n", 13);
    l.message(error, "Something failed!");
}