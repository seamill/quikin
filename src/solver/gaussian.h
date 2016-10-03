#ifndef _qk_solver_gaussian_H
#define _qk_solver_gaussian_H

// QK includes
#include "solver/solver.h"
#include "variable/variable_manager.h"

namespace qk
{

namespace solver
{

class gaussian: public solver
{
public:
    gaussian();
    ~gaussian();

    void setup(const double amplitude, const std::vector<double> & average, const double standard_deviation);

    void
    solve(const double time, qk::variable::variable_manager & variable_manager) const;

protected:

    double _amplitude;
    std::vector<double> _average;
    double _standard_deviation;

};

}
}

#endif // _qk_solver_gaussian_H
