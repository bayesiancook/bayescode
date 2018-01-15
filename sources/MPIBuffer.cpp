
#include "MPIBuffer.hpp"

template <>
MPIBuffer& operator<<(MPIBuffer& buffer, const double& t) {
    buffer.PutDouble(t);
    return buffer;
};

template <>
const MPIBuffer& operator>>(const MPIBuffer& buffer, double& t) {
    buffer.GetDouble(t);
    return buffer;
}

template <>
MPIBuffer& operator<<(MPIBuffer& buffer, const int& t) {
    buffer.PutInt(t);
    return buffer;
};

template <>
const MPIBuffer& operator>>(const MPIBuffer& buffer, int& t) {
    buffer.GetInt(t);
    return buffer;
}

template <>
unsigned int MPISize(const double& d) {
    return 1;
}

template <>
unsigned int MPISize(const int& i) {
    return 1;
}

template <>
double& operator+=(double& d, const MPIBuffer& buffer) {
    double c;
    buffer >> c;
    d += c;
    return d;
}

template <>
int& operator+=(int& i, const MPIBuffer& buffer) {
    int j;
    buffer >> j;
    i += j;
    return i;
}

template <>
vector<int>& operator+=(vector<int>& v, const MPIBuffer& buffer) {
    for (unsigned int i = 0; i < v.size(); i++) {
        v[i] += buffer;
    }
    return v;
}

template <>
vector<double>& operator+=(vector<double>& v, const MPIBuffer& buffer) {
    for (unsigned int i = 0; i < v.size(); i++) {
        v[i] += buffer;
    }
    return v;
}
