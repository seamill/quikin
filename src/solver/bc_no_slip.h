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

class bc_no_slip: public solver
{
public:
    bc_no_slip()=default;
    ~bc_no_slip()=default;

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
