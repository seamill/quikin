#ifndef _qk_solver_solver_H
#define _qk_solver_solver_H

// STL includes
#include <string>
#include <vector>

// QK includes
#include "variable/variable_id.h"
#include "variable/variable_manager.h"

namespace qk
{
namespace solver
{

class solver
{
public:

    solver()
    {

    }

    virtual ~solver()
    {

    }

    virtual void
    set_input_variables(const std::vector<qk::variable::variable_id> & input_variable_ids)
    {
        _input_variable_ids = input_variable_ids;
    }

    virtual void
    set_output_variables(const std::vector<qk::variable::variable_id> & output_variable_ids)
    {
        _output_variable_ids = output_variable_ids;
    }

    virtual void
    solve(qk::variable::variable_manager & variable_manager, const int tag) const = 0;

protected:

    std::vector<qk::variable::variable_id> _input_variable_ids;
    std::vector<qk::variable::variable_id> _output_variable_ids;

};

}
}

#endif // _qk_solver_solver_H
