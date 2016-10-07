#include "magnetic_pulse.h"

// STL include
#include <cmath>

// QK includes
#include "grid/grid.h"
#include "variable/variable.h"

namespace qk
{
namespace solver
{

void magnetic_pulse::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{

    const qk::grid::grid & grid = variable_manager.grid();

    double x[3];
    const double Bz = _amplitude;

    for (int i = 0; i < _output_variable_ids.size(); i++) {
        qk::variable::variable & var = variable_manager.output_variable(_output_variable_ids[i]);

        for (qk::indexer data_idx = var.indexer(); data_idx.exists(); data_idx.next()) {
            qk::data::extended_datachunk & data = var[data_idx];
            data.fill(0.);

            qk::range rho_range = data.internal_range();
            rho_range.set(data.num_dims()-1,5,6);

            for (qk::indexer idx = data.indexer(rho_range); idx.exists(); idx.next()) {
                grid.xc(idx, x);

                double r2 = 0.5*(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/_standard_deviation;

//                data[idx] = (r2<0.25) ? Bz : 0.;

                data[idx] = Bz*std::exp(-r2);

            }

        }

    }

}

}
}
