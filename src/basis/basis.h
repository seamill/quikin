#ifndef _qk_basis_basis_H
#define _qk_basis_basis_H

// STL includes
#include <string>

namespace qk
{
namespace basis
{

class basis
{
public:

    basis() :
            _name("void"),
            _num_points(0)
    {

    }

    basis(const std::string & name, int num_points) :
            _name(name),
            _num_points(num_points)
    {

    }

    virtual ~basis()
    {

    }

    int num_points() const
    {
        return _num_points;
    }

    const std::string & name() const
    {
        return _name;
    }

    friend bool operator <(const basis & a, const basis & b)
    {
        // Test if they have the same name
        if (a._name < b._name) {
            return true;
        }

        return a._num_points < b._num_points;

    }

    friend bool operator ==(const basis & a, const basis & b)
    {
        // Test if they have the same name
        if (a._name != b._name) {
            return false;
        }

        return a._num_points == b._num_points;
    }

    friend bool operator !=(const basis & a, const basis & b)
    {
        return !(a==b);
    }

protected:

    std::string _name;
    int _num_points;

};

}
}

#endif // _qk_basis_basis_H
