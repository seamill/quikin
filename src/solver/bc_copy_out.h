#ifndef _qk_solver_bc_copy_out_H
#define _qk_solver_bc_copy_out_H

// QK includes
#include "solver/solver.h"

// STL includes
#include <vector>

namespace qk
{

namespace solver
{

class bc_copy_out: public solver
{
public:
    bc_copy_out()=default;
    ~bc_copy_out()=default;

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

#endif // _qk_solver_bc_copy_out_H
