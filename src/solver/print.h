#ifndef _qk_solver_print_H
#define _qk_solver_print_H

// QK includes
#include "solver/solver.h"

namespace qk
{

namespace solver
{

class print: public solver
{
public:
    print() = default;
    ~print() = default;

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag = 0) const;

protected:

};

}
}

#endif // _qk_solver_print_H
