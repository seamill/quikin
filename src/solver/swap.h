#ifndef _qk_solver_swap_H
#define _qk_solver_swap_H

// QK includes
#include "solver/solver.h"
#include "variable/variable_manager.h"

namespace qk
{

namespace solver
{

class swap: public solver
{
public:
    swap();
    ~swap();

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag=0) const;

protected:

};

}
}

#endif // _qk_solver_swap_H
