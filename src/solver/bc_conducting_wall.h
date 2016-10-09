#ifndef _qk_solver_bc_conducting_wall_H
#define _qk_solver_bc_conducting_wall_H

// QK includes
#include "solver/solver.h"

// STL includes
#include <vector>

namespace qk
{

namespace solver
{

class bc_conducting_wall: public solver
{
public:
    bc_conducting_wall()=default;
    ~bc_conducting_wall()=default;

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

#endif // _qk_solver_bc_conducting_wall_H
