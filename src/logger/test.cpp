#include "Logger.hpp"

int main() {
    MessageFormat error(
        {bold("========================================================================\n"),
            bold_red("ERROR"), bold(" | ")},
        {bold("\n========================================================================\n")},
        {bold("      | ")});

    Logger l(std::cout);
    l.print("Hello %d!\n", 13);
    l.message(error, "Something failed!\nbadly!\nlike, very badly!");
}