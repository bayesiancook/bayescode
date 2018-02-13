
#ifndef NODEARRAY_H
#define NODEARRAY_H

#include "Tree.hpp"
#include <vector>
#include "MPIBuffer.hpp"

/**
 * \brief An interface for an abstract tree-structured array (indexed by branches) of constant references over objects of type T
 *
 * This abstract class is meant as a general read-only interface
 * returning a reference over a T object for any branch of a given phylogenetic tree (through the GetVal(int index) method).
 * Nodees are indexed by integers (ranging from 0 to GetNnode()-1).
 * This class, and its various implementations/specializations, follows the same logic as the classes deriving from Array,
 * except for the fact that the Node selectors are meant to define distributions across branches.
 * 
 */

template<class T> class NodeSelector	{

	public:
	virtual ~NodeSelector() {}

    //! return the number of branches of the underlying tree
	int GetNnode() const {return GetTree().GetNnode();}

    //! return a const reference to the underlying tree
	virtual const Tree& GetTree() const = 0;
    //! const access to array element by index
	virtual const T& GetVal(int index) const = 0;

    //! return size of entire array, when put into an MPI buffer
    unsigned int GetMPISize() const {return this->GetNnode() * MPISize(this->GetVal(0));}

    //! write array into MPI buffer
    void MPIPut(MPIBuffer& buffer) const {
        for (int i=0; i<this->GetNnode(); i++)  {
            buffer << this->GetVal(i);
        }
    }

    //! write array into generic output stream
    void ToStream(ostream& os) const    {
        for (int i=0; i<this->GetNnode(); i++)  {
            os << this->GetVal(i) << '\t';
        }
    }
};

/**
 * \brief An interface for an abstract tree-structured array (indexed by branches) of references over objects of type T
 *
 * Deriving from NodeSelector, NodeArray gives const access through GetVal(int index).
 * In addition, it provides a non-const access to through the [] operator.
 * In practice, this class is (should be?) used only for implementing arrays of all-distinct instances of objects of type T.
 * See also Array (and Selector) for a similar class hierarchy.
 */

template<class T> class NodeArray : public NodeSelector<T>	{

	public:
	virtual ~NodeArray() {}

    //! element-by-element copy (arrays should be of same size)
    void Copy(const NodeSelector<T>& from)  {

        if (this->GetNnode() != from.GetNnode())    {
            cerr << "error: branch arrays do not have same size\n";
            exit(1);
        }
        for (int i=0; i<this->GetNnode(); i++)  {
            (*this)[i] = from.GetVal(i);
        }
    }

    //! non-const access to array element by index
	virtual T& operator[](int index) = 0;

    //! get array from MPI buffer
    void MPIGet(const MPIBuffer& buffer)    {
        for (int i=0; i<this->GetNnode(); i++)  {
            buffer >> (*this)[i];
        }
    }

    //! get array from generic input stream
    void FromStream(istream& is)    {
        for (int i=0; i<this->GetNnode(); i++)  {
            is >> (*this)[i];
        }
    }
};

/**
 * \brief output stream operator for a NodeSelector<T>
 * 
 * in practice, calls the ToStream method of NodeSelector<T>
 */

template<class T> ostream& operator<<(ostream& os, const NodeSelector<T>& array)  {
    array.ToStream(os);
    return os;
}

/**
 * \brief input stream operator for a NodeArray<T>
 *
 * in practice, calls the FromStream method of NodeArray<T>
 */

template<class T> istream& operator>>(istream& is, NodeArray<T>& array) {
    array.FromStream(is);
    return is;
}

/**
 * \brief A NodeSelector<T> that returns a reference to the same value T for any branch
 *
 * Useful for implementing simple models assuming, e.g. the same substitution matrix for all branches.
 */

template<class T> class NodeHomogeneousSelector : public NodeSelector<T>	{

	public:
    //! \brief Constructor, taking as its arguments the tree and the value to be returned for any branch
	NodeHomogeneousSelector(const Tree* intree, const T& invalue) : tree(intree), value(invalue) {}
	~NodeHomogeneousSelector() {}

	const Tree& GetTree() const /*override*/ {return tree;}

    //! return a reference to the same value (i.e. value) for any index
	const T& GetVal(int index) const /*override*/ {return value;}

	private:
	const Tree& tree;
	const T& value;
};

/**
 * \brief The 'standard' implementation of an NodeArray<T>, simply as a std::vector<T>
 *
 */

template<class T> class SimpleNodeArray : public NodeArray<T>	{

	public:
    //! Constructor (with only the tree given as argument)
	SimpleNodeArray(const Tree& intree) : tree(intree), array(intree.GetNnode()) {}

    //! Constructor with tree and initializer value
	SimpleNodeArray(const Tree& intree, const T& initval) : tree(intree), array(intree.GetNnode(),initval) {}
	virtual ~SimpleNodeArray() {}

	const Tree& GetTree() const /*override*/ {return tree;}
	T& operator[](int index) /*override*/ {return array[index];}
	const T& GetVal(int index) const /*override*/ {return array[index];}

	protected:
	const Tree& tree;
	vector<T> array;
};

#endif
