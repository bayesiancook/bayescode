#include "Logger.hpp"

int main() {
    MessageFormat error(
        {bold_token(std::string(75, '=')), "\n", bold_red_token("ERROR"), bold_token(" | ")},
        {"\n", bold_token(std::string(75, '=')), "\n"}, {bold_token("      | ")});
    MessageFormat warning({Token({1, 32}, "WARNING"), bold_token(" ~ ")}, {bold_token(" ~\n")},
        {bold_token("         ")});
    MessageFormat info({Token({1, 36}, "INFO"), bold_token(" ("), Token(timestamp_date),
                           bold_token(", "), Token(timestamp_time), bold_token(") ~ ")},
        {bold_token(" ~\n")}, {Token(" ")});
    MessageFormat debug({Token({1, 35}, "DEBUG"), bold_token(" ~ ")}, {bold_token(" ~\n")},
        {bold_token("        ")});

    Logger l(std::cout);
    l.message(error, "Something failed!\nbadly!\nlike, very badly!");
    l.message(warning, "This is baaaad!\n baaad!");
    l.message(info, "This is a message %d", 3);
    std::cin.ignore();
    l.message(info, "This is a message");
    l.message(debug, "This is a\ndebug message");
}