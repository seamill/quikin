#ifndef _qk_solver_radial_shock_tube_H
#define _qk_solver_radial_shock_tube_H

// QK includes
#include "solver/solver.h"

namespace qk
{

namespace solver
{

class radial_shock_tube: public solver
{
public:
    radial_shock_tube();
    ~radial_shock_tube();

    void setup(const double r, const double density_jump, const double pressure_jump, const double gamma, const std::vector<double> & center = std::vector<double>())
    {
        _radius = r;
        _rho_jump = density_jump;
        _P_jump = pressure_jump;
        _gamma = gamma;
        _center = center;
        _center.resize(3,0.);
    }

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag = 0) const;

protected:

    std::vector<double> _center;
    double _radius;
    double _rho_jump;
    double _P_jump;
    double _gamma;

};

}
}

#endif // _qk_solver_radial_shock_tube_H
