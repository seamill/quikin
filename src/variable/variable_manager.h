#ifndef _qk_variable_variable_manager_H
#define _qk_variable_variable_manager_H

// STL includes
#include <string>
#include <vector>
#include <map>

// QK includes
#include "grid/grid.h"
#include "variable/variable_id.h"

namespace qk
{
namespace variable
{

class variable;

class variable_manager
{
public:

    variable_manager(const qk::grid::grid & grid);

    ~variable_manager();

    const qk::variable::variable &
    input_variable(const qk::variable::variable_id & var) const;

    qk::variable::variable &
    output_variable(const qk::variable::variable_id & var);

    const qk::grid::grid & grid() const
    {
        return _grid;
    }

    void
    write_VTK(const std::string & prefix, const std::string & suffix, const std::vector<qk::variable::variable_id> & variables) const;

protected:

    void
    add_variable(const qk::variable::variable_id & var);

    const qk::grid::grid & _grid;
    std::map<const qk::variable::variable_id, qk::variable::variable *> _variables;

};

}
}

#endif // _qk_variable_variable_manager_H
