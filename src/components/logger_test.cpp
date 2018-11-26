#include "Logger.hpp"

int main() {
    MessageFormat error({bold(std::string(75, '=')), "\n", bold_red("ERROR"), bold(" | ")},
        {"\n", bold(std::string(75, '=')), "\n"}, {bold("      | ")});
    MessageFormat warning(
        {Token({1, 32}, "WARNING"), bold(" ~ ")}, {bold(" ~\n")}, {bold("         ")});
    MessageFormat info({Token({1, 36}, "INFO"), bold(" ("), Token(timestamp), bold(") ")},
        {bold(" ~\n")}, {bold("       ")});
    MessageFormat debug({Token({1, 35}, "DEBUG"), bold(" ~ ")}, {bold(" ~\n")}, {bold("        ")});

    Logger l(std::cout);
    l.message(error, "Something failed!\nbadly!\nlike, very badly!");
    l.message(warning, "This is baaaad!\n baaad!");
    l.message(info, "This is a\nmessage");
    l.message(debug, "This is a\ndebug message");
}