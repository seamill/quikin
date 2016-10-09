#ifndef _qk_solver_rk1_H
#define _qk_solver_rk1_H

// QK includes
#include "solver/solver.h"

namespace qk
{

namespace solver
{

class rk1: public solver
{
public:

    enum
    {
        STAGE_0 = 0
    };

    rk1() = default;
    ~rk1() = default;

    int num_stages() const {return 1;}

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag = 0) const;

protected:

    void
    stage_0(const double dt,
        const qk::data::extended_datachunk & n0,
        const qk::data::extended_datachunk & rhs,
        qk::data::extended_datachunk & n1) const;

};

}
}

#endif // _qk_solver_rk1_H
