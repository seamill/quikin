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
    write_vtk(const std::string & prefix, const std::string & suffix, const std::vector<qk::variable::variable_id> & variables) const;

    void set_time(const double time) {_time=time;}
    void set_dt(const double dt) {_dt=dt;}

    double time() const {return _time;}
    double dt() const {return _dt;}

protected:

    void
    add_variable(const qk::variable::variable_id & var);

    const qk::grid::grid & _grid;
    std::map<const qk::variable::variable_id, qk::variable::variable *> _variables;

    double _dt;
    double _time;

};

}
}

#endif // _qk_variable_variable_manager_H
