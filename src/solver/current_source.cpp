#include "current_source.h"

// QK includes
#include "lib/exception.h"
#include "variable/variable_id.h"
#include "variable/variable.h"

// STL include
#include <cmath>
#include <iostream>
#include <omp.h>


namespace qk
{
namespace solver
{

namespace current_source_solve
{

void
fv_solve_idx(const double multiplier, const qk::data::extended_datachunk & fluid_data, qk::data::extended_datachunk & rhs_data)
{

    qk::range mom_rng = fluid_data.internal_range();
    mom_rng.set(mom_rng.num_dims()-1,1,4);

    qk::range E_rng = rhs_data.internal_range();
    E_rng.set(E_rng.num_dims()-1,0,3);

    qk::indexer p_idx = fluid_data.indexer(mom_rng);
    qk::indexer E_idx = rhs_data.indexer(E_rng);

    while(p_idx.exists() and E_idx.exists()){
        rhs_data[E_idx] += multiplier * fluid_data[p_idx];

        ++E_idx;
        ++p_idx;
    }

}

void
fv_solve(const double multiplier, const qk::data::extended_datachunk & fluid_data, qk::data::extended_datachunk & rhs_data)
{

    qk::range mom_rng = fluid_data.internal_range();
    mom_rng.set(mom_rng.num_dims()-1,0,1);

    qk::range E_rng = rhs_data.internal_range();
    E_rng.set(E_rng.num_dims()-1,0,1);

    qk::indexer p_idx = fluid_data.indexer(mom_rng);
    qk::indexer E_idx = rhs_data.indexer(E_rng);

    const int imin = p_idx.linear_index()/5;
    const int imax = (p_idx.final_linear_index()+1)/5;

    const double * __restrict__ p = fluid_data.data(p_idx) + 1;
    double * __restrict__ rhs = rhs_data.data(E_idx);
#pragma omp parallel for
    for(int i=imin;i<imax;++i){

        rhs[0] += multiplier*p[0];
        rhs[1] += multiplier*p[1];
        rhs[2] += multiplier*p[2];

        p+=5;
        rhs+=6;
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
        const qk::variable::variable & fluid = variable_manager.input_variable(_input_variable_ids[i]);
        for (qk::indexer idx = fluid.indexer(); idx.exists(); ++idx) {
            current_source_solve::fv_solve_idx(_multipliers[i],fluid[idx],rhs[idx]);
        }
    }

}

}
}
