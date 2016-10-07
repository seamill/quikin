#ifndef _qk_solver_ssprk3_H
#define _qk_solver_ssprk3_H

// QK includes
#include "solver/solver.h"
#include "grid/rectilinear.h"
#include "variable/variable_id.h"
#include "variable/variable_manager.h"

namespace qk
{

namespace solver
{

class ssprk3: public solver
{
public:

    enum
    {
        STAGE_0 = 0, STAGE_1 = 1, STAGE_2 = 2
    };

    ssprk3();
    ~ssprk3();

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag = 0) const;

protected:

    void
    stage_0(const double dt,
        const qk::data::extended_datachunk & n0,
        const qk::data::extended_datachunk & rhs,
        qk::data::extended_datachunk & n1) const;

    void
    stage_1(const double dt,
        const qk::data::extended_datachunk & n0,
        const qk::data::extended_datachunk & n1,
        const qk::data::extended_datachunk & rhs,
        qk::data::extended_datachunk & n2) const;

    void
    stage_2(const double dt,
        const qk::data::extended_datachunk & n0,
        const qk::data::extended_datachunk & n1,
        const qk::data::extended_datachunk & n2,
        const qk::data::extended_datachunk & rhs,
        qk::data::extended_datachunk & n3) const;

    qk::grid::rectilinear _grid;
    std::vector<double> _velocity;
};

}
}

#endif // _qk_solver_ssprk3_H
