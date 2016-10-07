#ifndef _qk_solver_current_source_H
#define _qk_solver_current_source_H

// QK includes
#include "solver/solver.h"

namespace qk
{

namespace solver
{

class current_source: public solver
{
public:

    current_source() = default;
    ~current_source() = default;

    void
    setup(const double eps0, const std::vector<double> & charges, const std::vector<double> & masses)
    {
        if(charges.size() != masses.size()){
            throw qk::exception("qk::solver::current_source::setup : Charges and masses are of different lengths.");
        }

        _multipliers.clear();
        for(int i=0;i<charges.size();i++){
            _multipliers.push_back( - charges[i]/masses[i]/eps0);
        }
    }

    void
    solve(qk::variable::variable_manager & variable_manager, const int tag=0) const;

protected:

    std::vector<double> _multipliers;

};

}
}

#endif // _qk_solver_current_source_H
