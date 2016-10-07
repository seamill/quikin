#ifndef _qk_solver_fill_H
#define _qk_solver_fill_H

// QK includes
#include "solver/solver.h"

namespace qk
{

namespace solver
{

class fill: public solver
{
public:
    fill() :
            _value(0)
    {

    }
    ~fill()
    {

    }

    void setup(const double value)
    {
        _value = value;
    }

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag = 0) const;

protected:

    double _value;

};

}
}

#endif // _qk_solver_gaussian_H
