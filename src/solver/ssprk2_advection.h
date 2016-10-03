#ifndef _qk_spatial_solver_advection_H
#define _qk_spatial_solver_advection_H

// QK includes
#include "solver/solver.h"
#include "grid/rectilinear.h"
#include "variable/variable_manager.h"

namespace qk
{

namespace solver
{

class ssprk2_advection: public solver
{
public:
    ssprk2_advection();
    ~ssprk2_advection();

    void setup(const std::vector<double> & velocity);

    void
    solve(const double time, qk::variable::variable_manager & variable_manager) const;

protected:

    qk::grid::rectilinear _grid;
    std::vector<double> _velocity;
};

}
}

#endif // _qk_spatial_solver_advection_H
