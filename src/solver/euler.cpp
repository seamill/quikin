#include "euler.h"

// QK includes
#include "lib/exception.h"
#include "grid/rectilinear.h"
#include "variable/variable_id.h"
#include "variable/variable.h"
#include "spatial_solver/functions.h"
#include "lib/parallelization.h"

// STL include
#include <cmath>
#include <iostream>

template<typename T>
T clamp(const T& n, const T& lower, const T& upper)
{
    return std::max(lower, std::min(n, upper));
}

namespace qk
{
namespace solver
{

namespace euler_solver
{

void solve_HLL_mp(const double gamma,
    const double dt,
    const qk::grid::rectilinear & grid,
    const qk::data::extended_datachunk & q_data,
    qk::data::extended_datachunk & rhs_data)
{

    // TODO : Make this scheme more efficient.
    // Due to issues with parallel write access, it is difficult to do this without calculating both left and right fluxes for each element
    // this can be countered by writing to a flux array and having a separate operation sum the fluxes properly into the rhs array
    // Right now we are double calculating all fluxes, which is a pretty big performance hit.

    double fl[5];
    double fr[5];
    double ql[5], qr[5], F[5];
    const double g1 = gamma - 1.;

    for (int dim = 0; dim < grid.num_dims(); ++dim) {

        const double d = grid.dx(dim);
        const int stride = q_data.stride(dim);
        qk::range itr_range(q_data.internal_range());
//        itr_range.expand(dim, -1, 0);
        itr_range.set(itr_range.num_dims() - 1, 0, 1);
        qk::indexer idx = q_data.indexer(itr_range);
        const int imin = idx.linear_index();
        const int imax = (idx.final_linear_index() + 1);

        const double * __restrict__ const q = q_data.data();
        double * __restrict__ const rhs = rhs_data.data();

#pragma omp parallel for private(fl,fr,ql,qr,F) num_threads(_QK_NUM_THREADS_)
        for (int i = imin; i < imax; i+=5) {

            const double * __restrict__ qc_ = q+i;
            const double * __restrict__ qn_ = qc_-stride;
            const double * __restrict__ qnn_ = qn_-stride;
            const double * __restrict__ qp_ = qc_+stride;
            const double * __restrict__ qpp_ = qp_+stride;
            double * __restrict__ rhs_ = rhs+i;

//          printf("thread %i : %i in %i, %i, stride %i\n",omp_get_thread_num(),i,imin,imax,stride);

            // First we do the left wall
            {
                // Reconstruct the conserved variables at the wall
                for (int j = 0; j < 5; ++j) {
                    const double & qn = qnn_[j];
                    const double & qc = qn_[j];
                    const double & qp = qc_[j];
                    const double & qpp = qp_[j];
                    const double rl = (qp - qc) / (qc - qn + 1.e-12);
                    const double rr = (qp - qc) / (qpp - qp + 1.e-12);
                    const double phil = (rl * rl + rl) / (rl * rl + 1);
                    const double phir = (rr * rr + rr) / (rr * rr + 1);
                    ql[j] = qc + 0.5 * phil * (qc - qn);
                    qr[j] = qp + 0.5 * phir * (qp - qpp);
                }

                // Calculate left and right pressures
                const double Pl = g1 * (ql[4] - 0.5 * (ql[1] * ql[1] + ql[2] * ql[2] + ql[3] * ql[3]) / ql[0]);
                const double Pr = g1 * (qr[4] - 0.5 * (qr[1] * qr[1] + qr[2] * qr[2] + qr[3] * qr[3]) / qr[0]);

                // Calculate speed of sound on each side of wall
                const double vsl = std::sqrt(gamma * Pl / ql[0]);
                const double vsr = std::sqrt(gamma * Pr / qr[0]);

                // Calculate face normal flow velocity
                const double unl = ql[dim + 1] / ql[0];
                const double unr = qr[dim + 1] / qr[0];

                // Calculate maximum left and right speeds
                const double sll = unl - vsl;
                const double srr = unr + vsr;

                // Calculate left flux
                fl[0] = ql[dim + 1];
                fl[1] = ql[1] * unl;
                fl[2] = ql[2] * unl;
                fl[3] = ql[3] * unl;
                fl[4] = (ql[4] + Pl) * unl;
                fl[dim + 1] += Pl;

                // Calculate right flux
                fr[0] = qr[dim + 1];
                fr[1] = qr[1] * unr;
                fr[2] = qr[2] * unr;
                fr[3] = qr[3] * unr;
                fr[4] = (qr[4] + Pr) * unr;
                fr[dim + 1] += Pr;

                // Create HLL flux multipliers - either 0 or 1
                const double fl_mult = double(sll > 0);
                const double fr_mult = double(srr < 0);
                const double fhll_mult = double(srr >= 0 and sll <= 0);

                // Use flux multipliers to define the actual flux
                for (int j = 0; j < 5; ++j) {
                    // This is the HLL flux for when the fluid is subsonic
                    F[j] = (srr * fl[j] - sll * fr[j] + sll * srr * (qr[j] - ql[j])) / (srr - sll);

                    // This is the full flux calculation with the flux multipliers
                    F[j] = fhll_mult * F[j] + fl_mult * fl[j] + fr_mult * fr[j];
                }

            }

            // We add the flux from the left wall
            for (int j = 0; j < 5; ++j) {
                rhs_[j] += F[j] / d;
            }

            // Now we do the right wall
            {
                // Reconstruct the conserved variables at the wall
                for (int j = 0; j < 5; ++j) {
                    const double & qn = qn_[j];
                    const double & qc = qc_[j];
                    const double & qp = qp_[j];
                    const double & qpp = qpp_[j];
                    const double rl = (qp - qc) / (qc - qn + 1.e-12);
                    const double rr = (qp - qc) / (qpp - qp + 1.e-12);
                    const double phil = (rl * rl + rl) / (rl * rl + 1);
                    const double phir = (rr * rr + rr) / (rr * rr + 1);
                    ql[j] = qc + 0.5 * phil * (qc - qn);
                    qr[j] = qp + 0.5 * phir * (qp - qpp);
                }

                // Calculate left and right pressures
                const double Pl = g1 * (ql[4] - 0.5 * (ql[1] * ql[1] + ql[2] * ql[2] + ql[3] * ql[3]) / ql[0]);
                const double Pr = g1 * (qr[4] - 0.5 * (qr[1] * qr[1] + qr[2] * qr[2] + qr[3] * qr[3]) / qr[0]);

                // Calculate speed of sound on each side of wall
                const double vsl = std::sqrt(gamma * Pl / ql[0]);
                const double vsr = std::sqrt(gamma * Pr / qr[0]);

                // Calculate face normal flow velocity
                const double unl = ql[dim + 1] / ql[0];
                const double unr = qr[dim + 1] / qr[0];

                // Calculate maximum left and right speeds
                const double sll = unl - vsl;
                const double srr = unr + vsr;

                // Calculate left flux
                fl[0] = ql[dim + 1];
                fl[1] = ql[1] * unl;
                fl[2] = ql[2] * unl;
                fl[3] = ql[3] * unl;
                fl[4] = (ql[4] + Pl) * unl;
                fl[dim + 1] += Pl;

                // Calculate right flux
                fr[0] = qr[dim + 1];
                fr[1] = qr[1] * unr;
                fr[2] = qr[2] * unr;
                fr[3] = qr[3] * unr;
                fr[4] = (qr[4] + Pr) * unr;
                fr[dim + 1] += Pr;

                // Create HLL flux multipliers - either 0 or 1
                const double fl_mult = double(sll > 0);
                const double fr_mult = double(srr < 0);
                const double fhll_mult = double(srr >= 0 and sll <= 0);

                // Use flux multipliers to define the actual flux
                for (int j = 0; j < 5; ++j) {
                    // This is the HLL flux for when the fluid is subsonic
                    F[j] = (srr * fl[j] - sll * fr[j] + sll * srr * (qr[j] - ql[j])) / (srr - sll);

                    // This is the full flux calculation with the flux multipliers
                    F[j] = fhll_mult * F[j] + fl_mult * fl[j] + fr_mult * fr[j];
                }

            }

            // We subtract the flux from the right wall
            for (int j = 0; j < 5; ++j) {
                rhs_[j] -= F[j] / d;
            }

        }
    }
}

void solve_HLL_fast(const double gamma,
    const double dt,
    const qk::grid::rectilinear & grid,
    const qk::data::extended_datachunk & q_data,
    qk::data::extended_datachunk & rhs_data)
{

    double fl[5];
    double fr[5];
    const double g1 = gamma - 1.;

    for (int dim = 0; dim < grid.num_dims(); dim++) {

        _QK_FV_UPWIND_RECON_SOLVE_START(dim, 5, q_data, rhs_data)
            const double Pl = g1 * (ql[4] - 0.5 * (ql[1] * ql[1] + ql[2] * ql[2] + ql[3] * ql[3]) / ql[0]);
            const double Pr = g1 * (qr[4] - 0.5 * (qr[1] * qr[1] + qr[2] * qr[2] + qr[3] * qr[3]) / qr[0]);

            const double vsl = std::sqrt(gamma * Pl / ql[0]);
            const double vsr = std::sqrt(gamma * Pr / qr[0]);

            const double unl = ql[dim + 1] / ql[0];
            const double unr = qr[dim + 1] / qr[0];

            const double sll = unl - vsl;
            const double srr = unr + vsr;

//            printf("(vsl,vsr)=(%1.2f,%1.2f), (sll,srr)=(%1.3f,%1.3f)\n",vsl,vsr,sll,srr);

            fl[0] = ql[dim + 1];
            fl[1] = ql[1] * unl;
            fl[2] = ql[2] * unl;
            fl[3] = ql[3] * unl;
            fl[4] = (ql[4] + Pl) * unl;
            fl[dim + 1] += Pl;

            fr[0] = qr[dim + 1];
            fr[1] = qr[1] * unr;
            fr[2] = qr[2] * unr;
            fr[3] = qr[3] * unr;
            fr[4] = (qr[4] + Pr) * unr;
            fr[dim + 1] += Pr;

            const double fl_mult = double(sll > 0);
            const double fr_mult = double(srr < 0);
            const double fhll_mult = double(srr >= 0 and sll <= 0);

            for (int j = 0; j < 5; ++j) {
                F[j] = (srr * fl[j] - sll * fr[j] + sll * srr * (qr[j] - ql[j])) / (srr - sll);

                F[j] = fhll_mult * F[j] + fl_mult * fl[j] + fr_mult * fr[j];
            }
        _QK_FV_UPWIND_RECON_SOLVE_END(5)
    }
}

}

void euler::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{
    if (_input_variable_ids.size() != 1) {
        throw qk::exception("qk::solver::euler::solve : Input must contain one variable (q).");
    }

    if (_output_variable_ids.size() != 1) {
        throw qk::exception("qk::solver::euler::solve : Output must contain one variable (rhs).");
    }

    const qk::grid::rectilinear & grid = dynamic_cast<const qk::grid::rectilinear &>(variable_manager.grid());

    const qk::variable::variable & q = variable_manager.input_variable(_input_variable_ids[0]);
    qk::variable::variable & rhs = variable_manager.output_variable(_output_variable_ids[0]);

    for (qk::indexer idx = q.indexer(); idx.exists(); ++idx) {
        euler_solver::solve_HLL_mp(_gamma, variable_manager.dt(), grid, q[idx], rhs[idx]);

    }

}

}
}
