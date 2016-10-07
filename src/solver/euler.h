#ifndef _qk_solver_euler_H
#define _qk_solver_euler_H

// QK includes
#include "solver/solver.h"

namespace qk
{

namespace solver
{

class euler: public solver
{
public:

    euler() = default;
    ~euler() = default;

    void
    setup(const double gamma)
    {
        _gamma=gamma;
    }

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag=0) const;

protected:

    double _gamma;

};

}
}

#endif // _qk_solver_euler_H
