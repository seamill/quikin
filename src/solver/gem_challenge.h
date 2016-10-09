#ifndef _qk_solver_gem_challenge_H
#define _qk_solver_gem_challenge_H

// QK includes
#include "solver/solver.h"

namespace qk
{

namespace solver
{

class gem_challenge: public solver
{
public:
    gem_challenge();
    ~gem_challenge();

    void setup(const double lambda,
        const double Lx,
        const double n0,
        const double n1,
        const double gamma,
        const double qe,
        const double me,
        const double ue,
        const double Te,
        const double qi,
        const double mi,
        const double ui,
        const double Ti,
        const double B0)
    {
        _lambda = lambda;
        _Lx = Lx;
        _n0 = n0;
        _n1 = n1;
        _gamma = gamma;
        _qe = qe;
        _qi = qi;
        _me = me;
        _mi = mi;
        _ue = ue;
        _ui = ui;
        _Te = Te;
        _Ti = Ti;
        _B0 = B0;
    }

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag = 0) const;

protected:

    double _lambda;
    double _Lx;
    double _n0;
    double _n1;
    double _gamma;
    double _qe;
    double _qi;
    double _me;
    double _mi;
    double _ue;
    double _ui;
    double _Te;
    double _Ti;
    double _B0;

};

}
}

#endif // _qk_solver_gem_challenge_H
