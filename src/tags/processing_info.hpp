#pragma once

#include <string>

template <class User, class Processing>
struct ProcessingInfo {
    using processing = Processing;

    User& user;
    std::string name;

    ProcessingInfo(User& user, std::string name = "") : user(user), name(name) {}
};

template <class Processing, class User>
auto make_processing_info(User& user) {
    return ProcessingInfo<User, Processing>(user);
}