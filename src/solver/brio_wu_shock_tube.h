#ifndef _qk_solver_brio_wu_shock_tube_H
#define _qk_solver_brio_wu_shock_tube_H

// QK includes
#include "solver/solver.h"

namespace qk
{

namespace solver
{

class brio_wu_shock_tube: public solver
{
public:
    brio_wu_shock_tube();
    ~brio_wu_shock_tube() = default;

    void setup(const double gamma,
        const double me,
        const double mi,
        const double nl,
        const double nr,
        const double Pl,
        const double Pr,
        const double Bx,
        const double By)
    {
        _gamma = gamma;
        _me = me;
        _mi = mi;
        _nl = nl;
        _nr = nr;
        _Pl = Pl;
        _Pr = Pr;
        _Bx = Bx;
        _By = By;

    }

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag = 0) const;

protected:

    double _gamma;
    double _me;
    double _mi;
    double _nl;
    double _nr;
    double _Pl;
    double _Pr;
    double _Bx;
    double _By;

};

}
}

#endif // _qk_solver_brio_wu_shock_tube_H
