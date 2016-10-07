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

void solve_rusanov(const double c,
    const qk::grid::rectilinear & grid,
    const qk::data::extended_datachunk & q_dat,
    qk::data::extended_datachunk & rhs_dat)
{

    const double * __restrict__ q = q_dat.data();
    double * __restrict__ rhs = rhs_dat.data();
    double ql[6];
    double qr[6];
    double F[6];

    const double absc = std::fabs(c);
    const double c2 = c * c;

    const double E_f_mult = 0.5 * c2;
    const double B_f_mult = 0.5;
    const double E_dq_mult = -0.5 * absc;
    const double B_dq_mult = -0.5 * absc;

    qk::range itr_range = q_dat.internal_range();
    itr_range.set(itr_range.num_dims() - 1, 0, 1);

    if (grid.num_dims() > 0) {

        const int s = q_dat.stride(0);
        const double dx = grid.dx(0);
        qk::range itr_range_2(itr_range);
        itr_range_2.expand(0, -1, 0);

//        std::cout << "Iteration " << itr_range_2 << std::endl;
        for (qk::indexer idx = q_dat.indexer(itr_range_2); idx.exists(); ++idx) {
            const int i = idx.linear_index();
            const int in = i - s;
            const int ip = i + s;
            const int ipp = ip + s;

            // Reconstruct surface values
            for (int j = 0; j < 6; ++j) {

                const double & qn = q[in + j];
                const double & qc = q[i + j];
                const double & qp = q[ip + j];
                const double & qpp = q[ipp + j];

                const double rl = (qp - qc) / (qc - qn + 1.e-12);
                const double rr = (qp - qc) / (qpp - qp + 1.e-12);

                const double phil = (rl * rl + rl) / (rl * rl + 1);
                const double phir = (rr * rr + rr) / (rr * rr + 1);

                ql[j] = qc + 0.5 * phil * (qc - qn);
                qr[j] = qp + 0.5 * phir * (qp - qpp);

            }

//            printf("index %i: \n\tql = [%1.3e, %1.3e, %1.3e, %1.3e, %1.3e, %1.3e] \n\tqr = [%1.3e, %1.3e, %1.3e, %1.3e, %1.3e, %1.3e]\n\tF = [%1.3e, %1.3e, %1.3e, %1.3e, %1.3e, %1.3e]\n",
//                i,
//                ql[0],ql[1],ql[2],ql[3],ql[4],ql[5],
//                qr[0],qr[1],qr[2],qr[3],qr[4],qr[5],
//                F[0],F[1],F[2],F[3],F[4],F[5]);

            F[0] = 0.;
            F[1] = -E_f_mult * (qr[5] + ql[5]) + E_dq_mult * (qr[1] - ql[1]);
            F[2] = E_f_mult * (qr[4] + ql[4]) + E_dq_mult * (qr[2] - ql[2]);
            ;
            F[3] = 0.;
            F[4] = B_f_mult * (qr[2] + ql[2]) + B_dq_mult * (qr[4] - ql[4]);
            F[5] = -B_f_mult * (qr[1] + ql[1]) + B_dq_mult * (qr[5] - ql[5]);

            for (int j = 0; j < 6; ++j) {
                rhs[i + j] -= F[j] / dx;
                rhs[ip + j] += F[j] / dx;
            }

        }

        if (grid.num_dims() > 1) {

            const int s = q_dat.stride(1);
            const double dx = grid.dx(1);
            qk::range itr_range_2(itr_range);
            itr_range_2.expand(1, -1, 0);

            //        std::cout << "Iteration " << itr_range_2 << std::endl;
            for (qk::indexer idx = q_dat.indexer(itr_range_2); idx.exists(); ++idx) {
                const int i = idx.linear_index();
                const int in = i - s;
                const int ip = i + s;
                const int ipp = ip + s;

                // Reconstruct surface values
                for (int j = 0; j < 6; ++j) {

                    const double & qn = q[in + j];
                    const double & qc = q[i + j];
                    const double & qp = q[ip + j];
                    const double & qpp = q[ipp + j];

                    const double rl = (qp - qc) / (qc - qn + 1.e-12);
                    const double rr = (qp - qc) / (qpp - qp + 1.e-12);

                    const double phil = (rl * rl + rl) / (rl * rl + 1);
                    const double phir = (rr * rr + rr) / (rr * rr + 1);

                    ql[j] = qc + 0.5 * phil * (qc - qn);
                    qr[j] = qp + 0.5 * phir * (qp - qpp);

                }

                //            printf("index %i: \n\tql = [%1.3e, %1.3e, %1.3e, %1.3e, %1.3e, %1.3e] \n\tqr = [%1.3e, %1.3e, %1.3e, %1.3e, %1.3e, %1.3e]\n\tF = [%1.3e, %1.3e, %1.3e, %1.3e, %1.3e, %1.3e]\n",
                //                i,
                //                ql[0],ql[1],ql[2],ql[3],ql[4],ql[5],
                //                qr[0],qr[1],qr[2],qr[3],qr[4],qr[5],
                //                F[0],F[1],F[2],F[3],F[4],F[5]);

                F[0] = E_f_mult * (qr[5] + ql[5]) + E_dq_mult * (qr[0] - ql[0]);
                F[1] = 0.;
                F[2] = -E_f_mult * (qr[3] + ql[3]) + E_dq_mult * (qr[2] - ql[2]);
                F[3] = -B_f_mult * (qr[2] + ql[2]) + B_dq_mult * (qr[3] - ql[3]);
                F[4] = 0.;
                F[5] = B_f_mult * (qr[0] + ql[0]) + B_dq_mult * (qr[5] - ql[5]);

                for (int j = 0; j < 6; ++j) {
                    rhs[i + j] -= F[j] / dx;
                    rhs[ip + j] += F[j] / dx;
                }

            }

        }

        if (grid.num_dims() > 2) {

            const int s = q_dat.stride(2);
            const double dx = grid.dx(2);
            qk::range itr_range_2(itr_range);
            itr_range_2.expand(2, -1, 0);

            //        std::cout << "Iteration " << itr_range_2 << std::endl;
            for (qk::indexer idx = q_dat.indexer(itr_range_2); idx.exists(); ++idx) {
                const int i = idx.linear_index();
                const int in = i - s;
                const int ip = i + s;
                const int ipp = ip + s;

                // Reconstruct surface values
                for (int j = 0; j < 6; ++j) {

                    const double & qn = q[in + j];
                    const double & qc = q[i + j];
                    const double & qp = q[ip + j];
                    const double & qpp = q[ipp + j];

                    const double rl = (qp - qc) / (qc - qn + 1.e-12);
                    const double rr = (qp - qc) / (qpp - qp + 1.e-12);

                    const double phil = (rl * rl + rl) / (rl * rl + 1);
                    const double phir = (rr * rr + rr) / (rr * rr + 1);

                    ql[j] = qc + 0.5 * phil * (qc - qn);
                    qr[j] = qp + 0.5 * phir * (qp - qpp);

                }

                //            printf("index %i: \n\tql = [%1.3e, %1.3e, %1.3e, %1.3e, %1.3e, %1.3e] \n\tqr = [%1.3e, %1.3e, %1.3e, %1.3e, %1.3e, %1.3e]\n\tF = [%1.3e, %1.3e, %1.3e, %1.3e, %1.3e, %1.3e]\n",
                //                i,
                //                ql[0],ql[1],ql[2],ql[3],ql[4],ql[5],
                //                qr[0],qr[1],qr[2],qr[3],qr[4],qr[5],
                //                F[0],F[1],F[2],F[3],F[4],F[5]);

                F[0] = -E_f_mult * (qr[4] + ql[4]) + E_dq_mult * (qr[0] - ql[0]);
                F[1] = E_f_mult * (qr[3] + ql[3]) + E_dq_mult * (qr[1] - ql[1]);
                F[2] = 0.;
                F[3] = B_f_mult * (qr[1] + ql[1]) + B_dq_mult * (qr[3] - ql[3]);
                F[4] = -B_f_mult * (qr[0] + ql[0]) + B_dq_mult * (qr[4] - ql[4]);
                F[5] = 0.;

                for (int j = 0; j < 6; ++j) {
                    rhs[i + j] -= F[j] / dx;
                    rhs[ip + j] += F[j] / dx;
                }

            }

        }

    }

}

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
            F[1] = -E_f_mult * (qr[5] + ql[5]) + E_dq_mult * (qr[1] - ql[1]);
            F[2] = E_f_mult * (qr[4] + ql[4]) + E_dq_mult * (qr[2] - ql[2]);
            F[3] = 0.;
            F[4] = B_f_mult * (qr[2] + ql[2]) + B_dq_mult * (qr[4] - ql[4]);
            F[5] = -B_f_mult * (qr[1] + ql[1]) + B_dq_mult * (qr[5] - ql[5]);
        _QK_FV_UPWIND_RECON_SOLVE_END(6)

    }

    if (grid.num_dims() > 1) {

        _QK_FV_UPWIND_RECON_SOLVE_START(1, 6, q_data, rhs_data)
            F[0] = E_f_mult * (qr[5] + ql[5]) + E_dq_mult * (qr[0] - ql[0]);
            F[1] = 0.;
            F[2] = -E_f_mult * (qr[3] + ql[3]) + E_dq_mult * (qr[2] - ql[2]);
            F[3] = -B_f_mult * (qr[2] + ql[2]) + B_dq_mult * (qr[3] - ql[3]);
            F[4] = 0.;
            F[5] = B_f_mult * (qr[0] + ql[0]) + B_dq_mult * (qr[5] - ql[5]);
        _QK_FV_UPWIND_RECON_SOLVE_END(6)
    }

    if (grid.num_dims() > 2) {

        _QK_FV_UPWIND_RECON_SOLVE_START(2, 6, q_data, rhs_data)
            F[0] = -E_f_mult * (qr[4] + ql[4]) + E_dq_mult * (qr[0] - ql[0]);
            F[1] = E_f_mult * (qr[3] + ql[3]) + E_dq_mult * (qr[1] - ql[1]);
            F[2] = 0.;
            F[3] = B_f_mult * (qr[1] + ql[1]) + B_dq_mult * (qr[3] - ql[3]);
            F[4] = -B_f_mult * (qr[0] + ql[0]) + B_dq_mult * (qr[4] - ql[4]);
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
