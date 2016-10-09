#include "gem_challenge.h"

// STL include
#include <cmath>

// QK includes
#include "grid/grid.h"
#include "variable/variable.h"

namespace qk
{
namespace solver
{

gem_challenge::gem_challenge()
{

}

gem_challenge::~gem_challenge()
{

}


void gem_challenge::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{

    double * __restrict__ e;
    double * __restrict__ i;
    double * __restrict__ f;
    double x[3] = {0};

    const double pi = 3.14159263;
    const double Lx = _Lx;
    const double Ly = 2*Lx;

    const double kx = pi/Lx;
    const double ky = 2*pi/Ly;

    const double ue = 2*_Te/_qe/_lambda/_B0;
    const double ui = 2*_Ti/_qi/_lambda/_B0;

    const double B1 = 0.1*_B0;

    qk::variable::variable & electrons = variable_manager.output_variable(_output_variable_ids[0]);
    qk::variable::variable & ions = variable_manager.output_variable(_output_variable_ids[1]);
    qk::variable::variable & field = variable_manager.output_variable(_output_variable_ids[2]);

    const qk::grid::grid & grid = variable_manager.grid();

    for (qk::indexer chunk_idx = electrons.indexer(); chunk_idx.exists(); chunk_idx.next()) {
        qk::data::extended_datachunk & e_chunk = electrons[chunk_idx];
        qk::data::extended_datachunk & i_chunk = ions[chunk_idx];
        qk::data::extended_datachunk & f_chunk = field[chunk_idx];
        e_chunk.fill(0.);
        i_chunk.fill(0.);
        f_chunk.fill(0.);

        qk::range e_rng = e_chunk.range();
        e_rng.set(e_rng.num_dims()-1,0,1);

        qk::range i_rng = i_chunk.range();
        i_rng.set(i_rng.num_dims()-1,0,1);

        qk::range f_rng = f_chunk.range();
        f_rng.set(f_rng.num_dims()-1,0,1);

        qk::indexer e_idx = e_chunk.indexer(e_rng);
        qk::indexer i_idx = i_chunk.indexer(i_rng);
        qk::indexer f_idx = f_chunk.indexer(f_rng);

        while(e_idx.exists() and i_idx.exists() and f_idx.exists()){
            e = e_chunk.data(e_idx);
            i = i_chunk.data(i_idx);
            f = f_chunk.data(f_idx);
            grid.xc(e_idx, x);

            const double sx = std::sin(kx*x[0]);
            const double cx = std::cos(kx*x[0]);
            const double sy = std::sin(ky*x[1]);
            const double cy = std::cos(ky*x[1]);
            double sech2 = 1. / std::cosh(x[0]/_lambda);
            sech2*=sech2;
            const double tanh = std::tanh(x[0]/_lambda);
            const double uvar = sech2 / (sech2 + _n1/_n0);

            const double n = _n0 * (sech2 + _n1/_n0);
            const double Pe = n * _Te;
            const double Pi = n * _Ti;

            e[0] = _me * n;
            i[0] = _mi * n;

            e[3] = e[0] * ue * uvar;
            i[3] = i[0] * ui * uvar;

            e[4] = Pe / (_gamma-1.) + 0.5 * e[3]*e[3]/e[0];
            i[4] = Pi / (_gamma-1.) + 0.5 * i[3]*i[3]/i[0];

            f[3] = B1 * ky*cx*sy;
            f[4] = _B0 * tanh - B1 * kx*sx*cy;
//            f[4] = - B1 * kx*sx*cy;

            ++e_idx;
            ++i_idx;
            ++f_idx;


        }
    }

}

}
}
