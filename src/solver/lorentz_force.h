#ifndef _qk_solver_lorentz_force_H
#define _qk_solver_lorentz_force_H

// QK includes
#include "solver/solver.h"

namespace qk
{

namespace solver
{

class lorentz_force: public solver
{
public:

    lorentz_force():
        _com(1.)
    {

    }
    ~lorentz_force() = default;

    void
    setup(const double charge, const double mass)
    {
        _com = charge / mass;
    }

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag=0) const;

protected:

    double _com;

};

}
}

#endif // _qk_solver_lorentz_force_H
