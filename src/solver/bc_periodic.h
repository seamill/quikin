#ifndef _qk_solver_bc_periodic_H
#define _qk_solver_bc_periodic_H

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
    bc_periodic();
    ~bc_periodic();

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag = 0) const;

protected:

};

}
}

#endif // _qk_solver_bc_periodic_H
