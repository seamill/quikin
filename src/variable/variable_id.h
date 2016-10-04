#ifndef _qk_variable_variable_id_H
#define _qk_variable_variable_id_H

// STL includes
#include <string>
#include <vector>

// QK includes
#include "basis/basis.h"
#include "data/dataset.h"
#include "lib/range.h"

namespace qk
{
namespace variable
{

class variable_id
{
public:

    variable_id();

    variable_id(const std::string & name, const std::vector<std::string> & component_names, const qk::basis::basis & basis);

    ~variable_id();

    variable_id unique_clone() const
    {
        variable_id vid(*this);

        std::stringstream ss;
        ss << "_clone_" << unique_tag();

        vid._name += ss.str();
        return vid;
    }

    static int unique_tag()
    {
        return _tag++;
    }

    const std::string & name() const {return _name;}
    const qk::basis::basis & basis() const {return _basis;}
    const std::vector<std::string> & component_names() const {return _component_names;}

    friend bool
    operator < (const variable_id & a, const variable_id & b)
    {
        // Test if they have the same name
        if(a._name < b._name){
            return true;
        } else if(a._name > b._name){
            return false;
        }

        if(a._component_names != b._component_names){
            throw qk::exception("qk::variable::variable_id::operator < : Shared variable name, but component mismatch.");
        }

        return a._basis < b._basis;

    }

    friend bool
    operator == (const variable_id & a, const variable_id & b)
    {
        // Test if they have the same name
        if(a._name != b._name){
            return false;
        }

        if(a._component_names != b._component_names){
            throw qk::exception("qk::variable::variable_id::operator < : Shared variable name, but component mismatch.");
        }

        return a._basis == b._basis;

    }

protected:

    static int _tag;
    std::string _name;
    qk::basis::basis _basis;
    std::vector<std::string> _component_names;

};

}
}

#endif // _qk_variable_variable_H
