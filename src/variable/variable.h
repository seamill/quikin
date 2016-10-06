#ifndef _qk_variable_variable_H
#define _qk_variable_variable_H

// STL includes
#include <string>
#include <vector>

// QK includes
#include "data/dataset.h"

namespace qk
{

namespace lib
{
class range;
}
namespace basis
{
class basis;
}

namespace variable
{

class variable: public qk::data::dataset
{
public:

    variable();

    variable(const std::string & name,
        const std::vector<std::string> & component_names,
        const qk::basis::basis & basis,
        const qk::range & global_mesh_range,
        const int num_ghost_layers=2);

    ~variable();

    void write_vtk(std::ofstream & file) const;

    const std::string & name() const
    {
        return _name;
    }

    const std::vector<std::string> & component_names() const
    {
        return _component_names;
    }

protected:

    std::string _name;
    std::vector<std::string> _component_names;

};

}
}

#endif // _qk_variable_variable_H
