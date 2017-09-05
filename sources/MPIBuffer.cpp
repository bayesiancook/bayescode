
#include "MPIBuffer.hpp"

template<> MPIBuffer& operator<<(MPIBuffer& buffer, const double& t)    {
    buffer.PutDouble(t);
    return buffer;
};

template<> const MPIBuffer& operator>>(const MPIBuffer& buffer, double& t)  {
    buffer.GetDouble(t);
    return buffer;
}

/*
template<> unsigned int MPISize<double>() {
    return 1;
};

template<> unsigned int MPISize<int>() {
    return 1;
};
*/

template<> unsigned int MPISize(const double& d)    {
    return 1;
}

template<> unsigned int MPISize(const int& i)   {
    return 1;
}

