#ifndef _qk_solver_magnetic_pulse_H
#define _qk_solver_magnetic_pulse_H

// QK includes
#include "solver/solver.h"

namespace qk
{

namespace solver
{

class magnetic_pulse: public solver
{
public:
    magnetic_pulse() = default;
    ~magnetic_pulse() = default;

    void setup(const double amplitude, const double standard_deviation)
    {
        _amplitude = amplitude;
        _standard_deviation = standard_deviation;
    }

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag = 0) const;

protected:

    double _amplitude;
    double _standard_deviation;

};

}
}

#endif // _qk_solver_magnetic_pulse_H
