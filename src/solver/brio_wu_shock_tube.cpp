#include "brio_wu_shock_tube.h"

// STL include
#include <cmath>

// QK includes
#include "grid/grid.h"
#include "variable/variable.h"

namespace qk
{
namespace solver
{

brio_wu_shock_tube::brio_wu_shock_tube() :
        _me(0.04),
        _mi(1.),
        _nl(1.),
        _nr(0.125),
        _Pl(1.),
        _Pr(0.1),
        _gamma(2.),
        _Bx(0.75),
        _By(1.)
{

}

void brio_wu_shock_tube::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{

    if (_output_variable_ids.size() != 3) {
        throw qk::exception("qk::solver::brio_wu_shock_tube::solve : Must have three output variables (electrons, ions, field).");
    }

    qk::variable::variable & electron_var = variable_manager.output_variable(_output_variable_ids[0]);
    qk::variable::variable & ion_var = variable_manager.output_variable(_output_variable_ids[1]);
    qk::variable::variable & field_var = variable_manager.output_variable(_output_variable_ids[2]);

    const qk::grid::grid & grid = variable_manager.grid();

    double x[3];
    const double P_mult = 1. / (_gamma - 1.);

    for (qk::indexer data_idx = electron_var.indexer(); data_idx.exists(); data_idx.next()) {
        qk::data::extended_datachunk & electrons = electron_var[data_idx];
        qk::data::extended_datachunk & ions = ion_var[data_idx];
        qk::data::extended_datachunk & field = field_var[data_idx];
        electrons.fill(0.);
        ions.fill(0.);
        field.fill(0.);

        qk::range rho_range = electrons.internal_range();
        rho_range.set(electrons.num_dims() - 1, 0, 1);

        for (qk::indexer idx = electrons.indexer(rho_range); idx.exists(); idx.next()) {
            grid.xc(idx, x);

            const double rho_e = _me * ((x[0] < 0.) ? _nl : _nr);
            const double rho_i = rho_e * _mi / _me;
            electrons[idx] = rho_e;
            ions[idx] = rho_i;
        }

        qk::range energy_range = electrons.internal_range();
        energy_range.set(electrons.num_dims() - 1, 4, 5);

        for (qk::indexer idx = electrons.indexer(energy_range); idx.exists(); idx.next()) {
            grid.xc(idx, x);

            const double P = ((x[0] < 0.) ? _Pl : _Pr) * P_mult;
            electrons[idx] = P;
            ions[idx] = P;
        }

        qk::range By_range = field.internal_range();
        By_range.set(electrons.num_dims() - 1, 4, 5);

        for (qk::indexer idx = field.indexer(By_range); idx.exists(); idx.next()) {
            grid.xc(idx, x);

            const double By = _By * ((x[0] < 0.) ? 1. : -1.) * P_mult;
            field[idx] = By;
        }

        qk::range Bx_range = field.internal_range();
        Bx_range.set(electrons.num_dims() - 1, 3, 4);
        field.fill(Bx_range, _Bx);

    }

}

}
}
