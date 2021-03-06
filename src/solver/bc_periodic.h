#ifndef _qk_solver_bc_periodic_H
#define _qk_solver_bc_periodic_H

// STL includes
#include <vector>

// QK includes
#include "solver/solver.h"
#include "variable/variable_manager.h"

namespace qk
{

namespace solver
{

class bc_periodic: public solver
{
public:
    bc_periodic()=default;
    ~bc_periodic()=default;

    void
    setup(const std::vector<int> & dims){
        _dims=dims;
    }

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag = 0) const;

protected:

    std::vector<int> _dims;

};

}
}

#endif // _qk_solver_bc_periodic_H
