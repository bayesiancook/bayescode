
#ifndef ARRAY_H
#define ARRAY_H

#include <vector>
#include "MPIBuffer.hpp"

/**
 * \brief An interface for an abstract indexed array of constant references over objects of type T
 *
 * This abstract class is meant as a general read-only interface
 * returning a reference over a T object for any index over a given range (through the GetVal(int index) method).
 * As an example, the PhyloProcess class requires a distribution of rates across sites, which is a Selector<double>.
 * In practice, the model could in fact assume one single rate for all sites (HomogeneousSelector<double>),
 * or a mixture of rates (MixtureSelector<double>), or an actual vector of rates (SimpleArray<double>),
 * all of which provide specific implementations of the generic interface provided by Selector<double>.
 * Note that all methods of Selector are const.
 * Thus, const Selector<T>& and Selector<T>& are equivalent in practice. However, const Selector<T>& is used everywhere
 * (to make it clear and explicit that no non-const access is intended, by the class receiving the Selector object).
 */

template<class T> class Selector	{

	public:
	virtual ~Selector() {}

    //! return array size
	virtual int GetSize() const = 0;
    //! const access to array element by index
	virtual const T& GetVal(int index) const = 0;

    //! return size of entire array, when put into an MPI buffer
    unsigned int GetMPISize() const {return this->GetSize() * MPISize(this->GetVal(0));}

    //! write array into MPI buffer
    void MPIPut(MPIBuffer& buffer) const {
        for (int i=0; i<this->GetSize(); i++)  {
            buffer << this->GetVal(i);
        }
    }

    //! write array into generic output stream
    void ToStream(ostream& os) const    {
        for (int i=0; i<this->GetSize(); i++)  {
            os << this->GetVal(i) << '\t';
        }
    }
};

/**
 * \brief An interface for an abstract indexed array of references over objects of type T
 *
 * Deriving from Selector, Array gives const access through GetVal(int index).
 * In addition, it provides a non-const access to through the [] operator.
 * In practice, this class is (should be?) used only for implementing arrays of all-distinct instances of objects of type T.
 * The SimpleArray class template provides a more explicit implementation (by encapsulating a std::vector<T>, see SimpleArray).
 * However, in some cases, it can be useful to implement an Array<T> as a std::vector<T*> instead, hence the specific 
 * declaration of this intermediate class.
 */

template<class T> class Array : public Selector<T>	{

	public:
	virtual ~Array() {}

    //! non-const access to array element by index
	virtual T& operator[](int index) = 0;

    //! element-by-element copy (arrays should be of same size)
    void Copy(const Selector<T>& from)  {

        if (this->GetSize() != from.GetSize())    {
            cerr << "error: arrays do not have same size\n";
            exit(1);
        }
        for (int i=0; i<this->GetSize(); i++)  {
            (*this)[i] = from.GetVal(i);
        }
    }

    //! get array from MPI buffer
    void MPIGet(const MPIBuffer& buffer)    {
        for (int i=0; i<this->GetSize(); i++)  {
            buffer >> (*this)[i];
        }
    }

    //! get array from generic input stream
    void FromStream(istream& is)    {
        for (int i=0; i<this->GetSize(); i++)  {
            is >> (*this)[i];
        }
    }
};

/**
 * \brief output stream operator for a Selector<T>
 * 
 * in practice, calls the ToStream method of Selector<T>
 */

template<class T> ostream& operator<<(ostream& os, const Selector<T>& array)  {
    array.ToStream(os);
    return os;
}

/**
 * \brief input stream operator for a Array<T>
 *
 * in practice, calls the FromStream method of Array<T>
 */

template<class T> istream& operator>>(istream& is, Array<T>& array) {
    array.FromStream(is);
    return is;
}

/**
 * \brief template for output stream operator for a std::vector<T>
 */

template<class T> ostream& operator<<(ostream& os, const vector<T>& array)  {
    for (unsigned int i=0; i<array.size(); i++) {
        os << array[i] << '\t';
    }
    return os;
}

/**
 * \brief template for input stream operator for a std::vector<T>
 */

template<class T> istream& operator>>(istream& is, vector<T>& array)    {
    for (unsigned int i=0; i<array.size(); i++) {
        is >> array[i];
    }
    return is;
}

/**
 * \brief A Selector<T> that returns a reference to the same value T for any index
 *
 * Useful for implementing simple models assuming, e.g. the same rate or the same substitution matrix for all sites.
 */

template<class T> class HomogeneousSelector : public Selector<T>	{

	public:
    //! \brief Constructor, taking as its arguments the size of the array and the value to be returned for any index
    //! \param insize: size of the array
    //! \param invalue: value to be returned for any index (type T should have a copy constructor)
	HomogeneousSelector(int insize, const T& invalue) : size(insize), value(invalue) {}

	~HomogeneousSelector() {}

	int GetSize() const override {return size;}

    //! return a reference to the same value (i.e. value) for any index
	const T& GetVal(int index) const override {return value;}

	private:
    //! array (abstract) size
	int size;

    //! value returned for any i=0..size-1
	const T& value;
};

/**
 * \brief The 'standard' implementation of an Array<T>, simply as a std::vector<T>
 *
 */

template<class T> class SimpleArray : public Array<T>	{

	public:
    //! Constructor (with only the array size given as argument)
	SimpleArray(int insize) : array(insize) {}

    //! Constructor with array size and initializer value
    SimpleArray(int insize, const T& initval) : array(insize,initval) {}

	virtual ~SimpleArray() {}

	int GetSize() const override {return array.size();}
	T& operator[](int index) override {return array[index];}
	const T& GetVal(int index) const override {return array[index];}

    //! return a const ref to the std::vector<T> of which this class is the interface
	const vector<T>& GetArray() const {return array;}

    //! Apply a permutation over the arguments, such that entry indexed by i after the call is equal to entry formely indexed by permut[i]
    virtual void Permute(const Selector<int>& permut) {
        if (permut.GetSize() != GetSize())  {
            cerr << "error in Array<T>::Permute: non matching array size\n";
            exit(1);
        }

        vector<T> tmp(GetSize(),GetVal(0));
        for (int i=0; i<GetSize(); i++) {
            tmp[i] = GetVal(permut.GetVal(i));
        }
        for (int i=0; i<GetSize(); i++) {
            (*this)[i] = tmp[i];
        }
    }

    //! Swap two entries
    void Swap(int cat1, int cat2)   {
        T tmp = (*this)[cat1];
        (*this)[cat1] = (*this)[cat2];
        (*this)[cat2] = tmp;
    }

	protected:
	vector<T> array;
};

/**
 * \brief A Selector that distributes the components of a mixture over an array of items, through a vector of allocations
 *
 * This Selector<T> maintains a pointer over an array of (say K) components and an array of (say N) integer allocation variables.
 * Then, for any index i, GetVal(i) returns a reference to the component to which item i is allocated (i.e. components[alloc[i]])
 */

template<class T> class MixtureSelector : public Selector<T>	{

	public:

    //! Constructor takes the array of components and the allocation vector
	MixtureSelector(const Selector<T>* incomponents, const Selector<int>* inalloc) : components(incomponents), alloc(inalloc)	{
	}
	~MixtureSelector() {}

    //! return size of array of components
	int GetSize() const override {return alloc->GetSize();}

    //! GetVal(i) (for i in 0..alloc->GetSize()-1) returns a reference to the component to which item i is allocated (i.e. components[alloc[i]])
	const T& GetVal(int i) const override {
		return components->GetVal(alloc->GetVal(i));
	}

	private:
	const Selector<T>* components;
	const Selector<int>* alloc;
};

#endif

