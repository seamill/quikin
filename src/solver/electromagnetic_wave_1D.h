#ifndef _qk_solver_electromagnetic_wave_1D_H
#define _qk_solver_electromagnetic_wave_1D_H

// QK includes
#include "solver/solver.h"

namespace qk
{

namespace solver
{

class electromagnetic_wave_1D: public solver
{
public:
    electromagnetic_wave_1D(){
        _dim=0;
        _k = 2*3.14159263;
        _c = 1.;
        _Bx=0;
        _By=1.;
        _Bz=1.;
    }
    ~electromagnetic_wave_1D() = default;

    void setup(const int dim, const double c, const double k, const double Bx, const double By, const double Bz)
    {
        _dim = dim;
        _c = c;
        _k = k;
        _Bx = Bx;
        _By = By;
        _Bz = Bz;
    }

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag = 0) const;

protected:

    int _dim;
    double _c;
    double _k;
    double _Bx, _By, _Bz;

};

}
}

#endif // _qk_solver_electromagnetic_wave_1D_H
