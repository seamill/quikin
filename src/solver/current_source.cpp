#include "current_source.h"

// STL include
#include <cmath>
#include <iostream>

// QK includes
#include "lib/exception.h"
#include "variable/variable_id.h"
#include "variable/variable.h"

namespace qk
{
namespace solver
{

namespace current_source_solve
{
void
fv_solve(const double multiplier, const qk::data::extended_datachunk & fluid_data, qk::data::extended_datachunk & rhs_data)
{

    qk::range int_rng = fluid_data.internal_range();
    int_rng.set(int_rng.num_dims()-1,0,1);
    qk::indexer idx = fluid_data.indexer(int_rng);
    const int imin = idx.linear_index();
    const int imax = idx.final_linear_index()+1;

    const double * __restrict__ fluid = fluid_data.data();
    double * __restrict__ rhs = rhs_data.data();

    const double * __restrict__ p;
    double * __restrict__ src;

    for(int i=imin;i<imax;++i){
        p = fluid + i*5 + 1;
        src = rhs + i*6;

        // We are solving
        // rhs Ex = (-q/e/eps0) * px
        // rhs Ey = (-q/e/eps0) * py
        // rhs Ez = (-q/e/eps0) * pz

        src[0] += multiplier * p[0];
        src[1] += multiplier * p[1];
        src[2] += multiplier * p[2];

    }
}

}


void current_source::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{
    if (_input_variable_ids.size() != _multipliers.size()) {
        throw qk::exception("qk::solver::current_source::solve : Input must contain same number of variables (fluids) as charges and masses.");
    }

    if (_output_variable_ids.size() != 1) {
        throw qk::exception("qk::solver::current_source::solve : Output must contain one variable (rhs).");
    }

    qk::variable::variable & rhs = variable_manager.output_variable(_output_variable_ids[0]);

    for(int i = 0; i < _multipliers.size(); ++i){
        const qk::variable::variable & fluid = variable_manager.input_variable(_input_variable_ids[0]);
        for (qk::indexer idx = fluid.indexer(); idx.exists(); ++idx) {
            current_source_solve::fv_solve(_multipliers[i],fluid[idx],rhs[idx]);
        }
    }

}

}
}
