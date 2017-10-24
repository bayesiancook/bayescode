
#ifndef MPIBUFFER_H
#define MPIBUFFER_H

#include <vector>
using namespace std;

class MPIBuffer {

    public:

    MPIBuffer(unsigned int insize) :buffer(0), size(insize), it(0) {
        buffer = new double[insize];
    }

    ~MPIBuffer() {
        delete[] buffer;
    }

    const double* GetBuffer() const {return buffer;}
    double* GetBuffer() {return buffer;}

    unsigned int GetSize() {return size;}

    template<class T> void Put(const T& t)  {
        t.MPIPut(*this);
    }

    template<class T> void Get(T& t) const {
        t.MPIGet(*this);
    }

    void PutDouble(const double& d)    {
        if (it == size) {
            // errror
        }
        buffer[it] = d;
        it++;
    }

    void GetDouble(double& d) const {
        if (it == size) {
            // errror
        }
        d = buffer[it];
        it++;
    }

    void PutInt(const int& i) {
        if (it == size) {
            // errror
        }
        buffer[it] = ((double) i);
        it++;
    }

    void GetInt(int& i) const {
        if (it == size) {
            // errror
        }
        double d = buffer[it];
        it++;
        i = (int) d;
    }

    private:

    double* buffer;
    unsigned int size;
    mutable unsigned int it;
};

// MPISize

template<class T> unsigned int MPISize(const T& t)  {
    return t.GetMPISize();
}

template<> unsigned int MPISize(const double& d);
template<> unsigned int MPISize(const int& i);

template<class T> unsigned int MPISize(const vector<T>& t)  {
    return t.size() * MPISize(t[0]);
}

// MPI operator << and >>

template<class T> MPIBuffer& operator<<(MPIBuffer& buffer, const T& t) {
    buffer.Put(t);
    return buffer;
};

template<class T> const MPIBuffer& operator>>(const MPIBuffer& buffer, T& t)    {
    buffer.Get(t);
    return buffer;
}

template<class T> MPIBuffer& operator<<(MPIBuffer& buffer, const vector<T>& v)  {
    for (unsigned int i=0; i<v.size(); i++) {
        buffer << v[i];
    }
    return buffer;
}

template<class T> const MPIBuffer& operator>>(const MPIBuffer& buffer, vector<T>& v)    {
    for (unsigned int i=0; i<v.size(); i++) {
        buffer >> v[i];
    }
    return buffer;
}

template<> MPIBuffer& operator<<(MPIBuffer& buffer, const double& t);
template<> const MPIBuffer& operator>>(const MPIBuffer& buffer, double& t);
template<> MPIBuffer& operator<<(MPIBuffer& buffer, const int& t);
template<> const MPIBuffer& operator>>(const MPIBuffer& buffer, int& t);

// MPI operator += with buffer

template<class T> T& operator +=(T& t, const MPIBuffer& buffer) {
    t.Add(buffer);
    return t;
}

template<> double& operator +=(double& d, const MPIBuffer& buffer);
template<> int& operator +=(int& i, const MPIBuffer& buffer);

template<> vector<int>& operator += (vector<int>& v, const MPIBuffer& buffer);
template<> vector<double>& operator += (vector<double>& v, const MPIBuffer& buffer);

#endif

