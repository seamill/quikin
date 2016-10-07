#include "radial_shock_tube.h"

// STL include
#include <cmath>

// QK includes
#include "grid/grid.h"
#include "variable/variable.h"

namespace qk
{
namespace solver
{

radial_shock_tube::radial_shock_tube():
    _radius(0.125),
    _rho_jump(0.125),
    _P_jump(0.1),
    _gamma(1.4)
{

}

radial_shock_tube::~radial_shock_tube()
{

}


void radial_shock_tube::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{

    const qk::grid::grid & grid = variable_manager.grid();

    double x[3];
    const double R2 = _radius*_radius;
    const double rho0 = 1.;
    const double e0 = 1. / (_gamma-1.);

    for (int i = 0; i < _output_variable_ids.size(); i++) {
        qk::variable::variable & var = variable_manager.output_variable(_output_variable_ids[i]);

        for (qk::indexer data_idx = var.indexer(); data_idx.exists(); data_idx.next()) {
            qk::data::extended_datachunk & data = var[data_idx];
            data.fill(0.);

            qk::range rho_range = data.internal_range();
            rho_range.set(data.num_dims()-1,0,1);

            for (qk::indexer idx = data.indexer(rho_range); idx.exists(); idx.next()) {
                grid.xc(idx, x);

                double r2 = (x[0]-_center[0])*(x[0]-_center[0]);
                r2 += (x[1]-_center[1])*(x[1]-_center[1]);
                r2 += (x[2]-_center[2])*(x[1]-_center[2]);

                const double rho = (r2 < R2) ? rho0 : rho0 * _rho_jump;
                data[idx] = rho;
            }

            qk::range e_range = data.internal_range();
            e_range.set(data.num_dims()-1,4,5);

            for (qk::indexer idx = data.indexer(e_range); idx.exists(); idx.next()) {
                grid.xc(idx, x);

                double r2 = (x[0]-_center[0])*(x[0]-_center[0]);
                r2 += (x[1]-_center[1])*(x[1]-_center[1]);
                r2 += (x[2]-_center[2])*(x[1]-_center[2]);

                const double e = (r2 < R2) ? e0 : e0 * _P_jump;
                data[idx] = e;
            }

        }

    }

}

}
}
