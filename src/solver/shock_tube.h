#ifndef _qk_solver_shock_tube_H
#define _qk_solver_shock_tube_H

// QK includes
#include "solver/solver.h"

namespace qk
{

namespace solver
{

class shock_tube: public solver
{
public:
    shock_tube();
    ~shock_tube();

    void setup(const double density_jump, const double pressure_jump, const double gamma)
    {
        _rho_jump = density_jump;
        _P_jump = pressure_jump;
        _gamma = gamma;
    }

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag = 0) const;

protected:

    double _rho_jump;
    double _P_jump;
    double _gamma;

};

}
}

#endif // _qk_solver_shock_tube_H
