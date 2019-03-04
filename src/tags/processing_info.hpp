#pragma once

template <class User, class Processing>
struct ProcessingInfo {
    using processing = Processing;

    User& user;

    ProcessingInfo(User& user) : user(user) {}
};

template <class Processing, class User>
auto make_processing_info(User& user) {
    return ProcessingInfo<User, Processing>(user);
}