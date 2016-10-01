#ifndef _qk_lib_exception_H
#define _qk_lib_exception_H

// STL includes
#include <exception>
#include <sstream>

namespace qk
{

class exception: public std::exception
{
public:

    exception()
    {

    }

    exception(const std::string error)
    {
        _what << error;
    }

    exception(const std::exception & other)
    {
        _what << std::string(other.what());
    }

    exception(const qk::exception & other)
    {
        _what << other._what.str();
    }

    virtual ~exception() throw()
    {

    }


    const char* what() const noexcept
    {
        return _what.str().c_str();
    }

    exception& operator=(const qk::exception & other)
    {
        this->_what.clear();
        *this << std::string(other.what());
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const qk::exception& qke)
    {
        os << qke._what.str();
        return os;
    }

    exception &
    operator<<(const std::string & s)
    {
        _what << s;
        return *this;
    }

protected:

    std::stringstream _what;
};

}

#endif // _qk_lib_exception_H

