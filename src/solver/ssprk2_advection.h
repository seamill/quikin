#ifndef _qk_spatial_solver_advection_H
#define _qk_spatial_solver_advection_H

// QK includes
#include "solver/solver.h"
#include "grid/rectilinear.h"
#include "variable/variable_id.h"
#include "variable/variable_manager.h"

namespace qk
{

namespace solver
{

class ssprk2_advection: public solver
{
public:

    enum{
        STAGE_0 = 0,
        STAGE_1 = 1
    };

    ssprk2_advection();
    ~ssprk2_advection();

    void
    setup(const std::vector<double> & velocity);

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag=0) const;

protected:

    void
    stage_0(const double dt, const qk::grid::rectilinear & grid, const qk::data::extended_datachunk & n0, qk::data::extended_datachunk & n1) const;

    void
    stage_1(const double dt, const qk::grid::rectilinear & grid, const qk::data::extended_datachunk & n0, const qk::data::extended_datachunk & n1, qk::data::extended_datachunk & n2) const;


    qk::grid::rectilinear _grid;
    std::vector<double> _velocity;
};

}
}

#endif // _qk_spatial_solver_advection_H
