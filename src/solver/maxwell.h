#ifndef _qk_solver_maxwell_H
#define _qk_solver_maxwell_H

// QK includes
#include "solver/solver.h"

namespace qk
{

namespace solver
{

class maxwell: public solver
{
public:

    maxwell() :
            _c(1.)
    {

    }
    ~maxwell() = default;

    void setup(const double c)
    {
        _c = c;
    }

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag = 0) const;

protected:

    double _c;

};

}
}

#endif // _qk_solver_maxwell_H
