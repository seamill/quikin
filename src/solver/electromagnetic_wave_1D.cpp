#include "electromagnetic_wave_1D.h"

// QK includes
#include "grid/grid.h"
#include "variable/variable.h"

// STL include
#include <cmath>

namespace qk
{
namespace solver
{

void electromagnetic_wave_1D::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{

    const qk::grid::grid & grid = variable_manager.grid();

    double x[3];
    double * __restrict__ f;

    for (const qk::variable::variable_id var_id : _output_variable_ids) {
        qk::variable::variable & var = variable_manager.output_variable(var_id);

        for (qk::indexer chunk_idx = var.indexer(); chunk_idx.exists(); chunk_idx.next()) {
            qk::data::extended_datachunk & chunk = var[chunk_idx];
            chunk.fill(0.);

            qk::range field_rng = chunk.range();
            field_rng.set(field_rng.num_dims()-1,0,1);

            if(_dim == 0){
                const double Ey = _c * _Bz;
                const double Ez = - _c * _By;

                for(qk::indexer idx = chunk.indexer(field_rng); idx.exists(); ++idx){
                    f = chunk.data(idx);
                    grid.xc(idx,x);
                    const double s = std::sin(_k * x[0]);

                    f[0] = 0;
                    f[1] = Ey * s;
                    f[2] = Ez * s;

                    f[3] = _Bx;
                    f[4] = _By * s;
                    f[5] = _Bz * s;

                }
            }

            if(_dim == 1){
                const double Ex = - _c * _Bz;
                const double Ez = _c * _Bx;

                for(qk::indexer idx = chunk.indexer(field_rng); idx.exists(); ++idx){
                    f = chunk.data(idx);
                    grid.xc(idx,x);
                    const double s = std::sin(_k * x[1]);

                    f[0] = Ex * s;
                    f[1] = 0.;
                    f[2] = Ez * s;

                    f[3] = _Bx * s;
                    f[4] = _By;
                    f[5] = _Bz * s;

                }
            }

            if(_dim == 2){
                const double Ex = _c * _By;
                const double Ey = - _c * _Bx;

                for(qk::indexer idx = chunk.indexer(field_rng); idx.exists(); ++idx){
                    f = chunk.data(idx);
                    grid.xc(idx,x);
                    const double s = std::sin(_k * x[2]);

                    f[0] = Ex * s;
                    f[1] = Ey * s;
                    f[2] = 0.;

                    f[3] = _Bx * s;
                    f[4] = _By * s;
                    f[5] = _Bz;

                }
            }

        }
    }

}

}
}
