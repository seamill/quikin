#include "lorentz_force.h"

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

namespace lorentz_force_solve
{
void
low_order_solve(const double com, const qk::data::extended_datachunk & fluid_data, const qk::data::extended_datachunk & field_data, qk::data::extended_datachunk & rhs_data)
{

    qk::range int_rng = fluid_data.internal_range();
    int_rng.set(int_rng.num_dims()-1,0,1);
    qk::indexer idx = fluid_data.indexer(int_rng);
    const int imin = idx.linear_index();
    const int imax = idx.final_linear_index()+1;

    const double * __restrict__ fluid = fluid_data.data();
    const double * __restrict__ field = field_data.data();
    double * __restrict__ rhs = rhs_data.data();

    const double * __restrict__ q;
    const double * __restrict__ p;
    const double * __restrict__ E;
    const double * __restrict__ B;
    double * __restrict__ r;

    for(int i=imin;i<imax;++i){
        q = fluid + i*5;
        p = q + 1;
        E = fluid + i*6;
        B = E+3;
        r = rhs + i*5;

        // We are solving
        // rhs px = com * ( rho * Ex + py * Bz - pz * By)
        // rhs py = com * ( rho * Ey + pz * Bx - px * Bz)
        // rhs pz = com * ( rho * Ez + px * By - py * Bx)
        // rhs e = com * ( px * Ex + py * Ey + pz * Ez)

        r[1] += com * (q[0]*E[0] + p[1]*B[2]-p[2]*B[1]);
        r[2] += com * (q[0]*E[1] + p[2]*B[0]-p[0]*B[2]);
        r[3] += com * (q[0]*E[2] + p[0]*B[1]-p[1]*B[0]);
        r[4] += com * (p[0]*E[0]+p[1]*E[1]+p[2]*E[2]);

    }
}

}


void lorentz_force::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{
    if (_input_variable_ids.size() != 2) {
        throw qk::exception("qk::solver::lorentz_force::solve : Input must contain two variables (fluid, field).");
    }

    if (_output_variable_ids.size() != 1) {
        throw qk::exception("qk::solver::lorentz_force::solve : Output must contain one variable (rhs).");
    }

    const qk::variable::variable & fluid = variable_manager.input_variable(_input_variable_ids[0]);
    const qk::variable::variable & field = variable_manager.input_variable(_input_variable_ids[1]);
    qk::variable::variable & rhs = variable_manager.output_variable(_output_variable_ids[0]);

    for (qk::indexer idx = fluid.indexer(); idx.exists(); ++idx) {
        lorentz_force_solve::low_order_solve(_com,fluid[idx],field[idx],rhs[idx]);

    }

}

}
}
