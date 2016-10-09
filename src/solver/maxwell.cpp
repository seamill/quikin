#include "maxwell.h"

// STL include
#include <cmath>

// QK includes
#include "lib/exception.h"
#include "grid/rectilinear.h"
#include "variable/variable_id.h"
#include "variable/variable.h"
#include "spatial_solver/functions.h"

namespace qk
{
namespace solver
{

namespace maxwell_solver
{

void solve_rusanov_streamlined(const double c,
    const qk::grid::rectilinear & grid,
    const qk::data::extended_datachunk & q_data,
    qk::data::extended_datachunk & rhs_data)
{

    const double absc = std::fabs(c);
    const double c2 = c * c;

    const double E_f_mult = 0.5 * c2;
    const double B_f_mult = 0.5;
    const double E_dq_mult = -0.5 * absc;
    const double B_dq_mult = -0.5 * absc;


    if (grid.num_dims() > 0) {

        _QK_FV_UPWIND_RECON_SOLVE_START(0, 6, q_data, rhs_data)
            F[0] = 0.;
            F[1] = E_f_mult * (qr[5] + ql[5]) + E_dq_mult * (qr[1] - ql[1]);
            F[2] = -E_f_mult * (qr[4] + ql[4]) + E_dq_mult * (qr[2] - ql[2]);
            F[3] = 0.;
            F[4] = -B_f_mult * (qr[2] + ql[2]) + B_dq_mult * (qr[4] - ql[4]);
            F[5] = B_f_mult * (qr[1] + ql[1]) + B_dq_mult * (qr[5] - ql[5]);
        _QK_FV_UPWIND_RECON_SOLVE_END(6)

    }

    if (grid.num_dims() > 1) {

        _QK_FV_UPWIND_RECON_SOLVE_START(1, 6, q_data, rhs_data)
            F[0] = -E_f_mult * (qr[5] + ql[5]) + E_dq_mult * (qr[0] - ql[0]);
            F[1] = 0.;
            F[2] = E_f_mult * (qr[3] + ql[3]) + E_dq_mult * (qr[2] - ql[2]);
            F[3] = B_f_mult * (qr[2] + ql[2]) + B_dq_mult * (qr[3] - ql[3]);
            F[4] = 0.;
            F[5] = -B_f_mult * (qr[0] + ql[0]) + B_dq_mult * (qr[5] - ql[5]);
        _QK_FV_UPWIND_RECON_SOLVE_END(6)
    }

    if (grid.num_dims() > 2) {

        _QK_FV_UPWIND_RECON_SOLVE_START(2, 6, q_data, rhs_data)
            F[0] = E_f_mult * (qr[4] + ql[4]) + E_dq_mult * (qr[0] - ql[0]);
            F[1] = -E_f_mult * (qr[3] + ql[3]) + E_dq_mult * (qr[1] - ql[1]);
            F[2] = 0.;
            F[3] = -B_f_mult * (qr[1] + ql[1]) + B_dq_mult * (qr[3] - ql[3]);
            F[4] = B_f_mult * (qr[0] + ql[0]) + B_dq_mult * (qr[4] - ql[4]);
            F[5] = 0.;
        _QK_FV_UPWIND_RECON_SOLVE_END(6)
    }

}

}

void maxwell::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{
    if (_input_variable_ids.size() != 1) {
        throw qk::exception("qk::solver::maxwell::solve : Input must contain one variable (q).");
    }

    if (_output_variable_ids.size() != 1) {
        throw qk::exception("qk::solver::maxwell::solve : Output must contain one variable (rhs).");
    }

    const qk::grid::rectilinear & grid = dynamic_cast<const qk::grid::rectilinear &>(variable_manager.grid());

    const qk::variable::variable & q = variable_manager.input_variable(_input_variable_ids[0]);
    qk::variable::variable & rhs = variable_manager.output_variable(_output_variable_ids[0]);

    for (qk::indexer idx = q.indexer(); idx.exists(); ++idx) {
        maxwell_solver::solve_rusanov_streamlined(_c, grid, q[idx], rhs[idx]);
    }

}

}
}
