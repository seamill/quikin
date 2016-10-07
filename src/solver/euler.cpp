#include "euler.h"

// STL include
#include <cmath>
#include <iostream>

// QK includes
#include "lib/exception.h"
#include "grid/rectilinear.h"
#include "variable/variable_id.h"
#include "variable/variable.h"
#include "spatial_solver/functions.h"

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

void solve(const double gamma,
    const qk::grid::rectilinear & grid,
    const qk::data::extended_datachunk & q_dat,
    qk::data::extended_datachunk & rhs_dat)
{

    const double * __restrict__ q = q_dat.data();
    double * __restrict__ rhs = rhs_dat.data();
    double qf[5];
    double flux[5];
    const double * p = qf + 1;
    const double g1 = gamma - 1.;

    qk::range itr_range = q_dat.internal_range();
    itr_range.set(itr_range.num_dims() - 1, 0, 1);

    for (int dim = 0; dim < grid.num_dims(); dim++) {
        const double dx = grid.dx(dim);
        qk::range itr_range_2(itr_range);
        itr_range_2.expand(dim, -1, 0);
        for (qk::indexer idx = q_dat.indexer(itr_range_2); idx.exists(); ++idx) {
            const int s = q_dat.stride(dim);
            const int i = idx.linear_index();
            const int in = i - s;
            const int ip = i + s;
            const int ipp = ip + s;
            // Reconstruct surface values
            for (int j = 0; j < 5; j++) {
                const double & qn = q[in + j];
                const double & qc = q[i + j];
                const double & qp = q[ip + j];
                const double & qpp = q[ipp + j];

                qf[j] = (7. / 12.) * (qc + qp) - (1. / 12.) * (qn + qpp);
            }

            const double & rho = qf[0];
            const double & e = qf[4];

            const double P = g1 * (e - 0.5 * (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) / rho);

            const double un = p[dim] / rho;

            flux[0] = p[dim];
            flux[1] = p[0] * un;
            flux[2] = p[1] * un;
            flux[3] = p[2] * un;
            flux[4] = (e + P) * un;
            flux[dim] += P;

            for (int j = 0; j < 5; ++j) {
                rhs[i + j] -= flux[j] / dx;
                rhs[ip + j] += flux[j] / dx;
            }

        }
    }
}

void solve_HLL(const double gamma,
    const double dt,
    const qk::grid::rectilinear & grid,
    const qk::data::extended_datachunk & q_dat,
    qk::data::extended_datachunk & rhs_dat)
{

    const double * __restrict__ q = q_dat.data();
    double * __restrict__ rhs = rhs_dat.data();
    double ql[5];
    double qr[5];
    double fl[5];
    double fr[5];
    const double * pl = ql + 1;
    const double * pr = qr + 1;
    const double g1 = gamma - 1.;

    qk::range itr_range = q_dat.internal_range();
    itr_range.set(itr_range.num_dims() - 1, 0, 1);

    for (int dim = 0; dim < grid.num_dims(); dim++) {
        const int s = q_dat.stride(dim);
        const double dx = grid.dx(dim);
        qk::range itr_range_2(itr_range);
        itr_range_2.expand(dim, -1, 0);
        for (qk::indexer idx = q_dat.indexer(itr_range_2); idx.exists(); ++idx) {
            const int i = idx.linear_index();
            const int in = i - s;
            const int ip = i + s;
            const int ipp = ip + s;

            // Reconstruct surface values
            for (int j = 0; j < 5; ++j) {

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

            const double Pl = g1 * (ql[4] - 0.5 * (pl[0] * pl[0] + pl[1] * pl[1] + pl[2] * pl[2]) / ql[0]);
            const double Pr = g1 * (qr[4] - 0.5 * (pr[0] * pr[0] + pr[1] * pr[1] + pr[2] * pr[2]) / qr[0]);

            const double vsl = std::sqrt(gamma * Pl / ql[0]);
            const double vsr = std::sqrt(gamma * Pr / qr[0]);

            const double unl = pl[dim] / ql[0];
            const double unr = pr[dim] / qr[0];

            const double sll = unl - vsl;
            const double srr = unr + vsr;

//            printf("(vsl,vsr)=(%1.2f,%1.2f), (sll,srr)=(%1.3f,%1.3f)\n",vsl,vsr,sll,srr);

            fl[0] = pl[dim];
            fl[1] = pl[0] * unl;
            fl[2] = pl[1] * unl;
            fl[3] = pl[2] * unl;
            fl[4] = (ql[4] + Pl) * unl;
            fl[dim + 1] += Pl;

            fr[0] = pr[dim];
            fr[1] = pr[0] * unr;
            fr[2] = pr[1] * unr;
            fr[3] = pr[2] * unr;
            fr[4] = (qr[4] + Pr) * unr;
            fr[dim + 1] += Pr;

            for (int j = 0; j < 5; ++j) {
                double F = 0;
                if (sll > 0) {
                    F = fl[j];
                } else if (srr < 0) {
                    F = fr[j];
                } else {
                    F = (srr * fl[j] - sll * fr[j] + sll * srr * (qr[j] - ql[j])) / (srr - sll);
                }

                rhs[i + j] -= F / dx;
                rhs[ip + j] += F / dx;
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

            const double unl = ql[dim+1] / ql[0];
            const double unr = qr[dim+1] / qr[0];

            const double sll = unl - vsl;
            const double srr = unr + vsr;

//            printf("(vsl,vsr)=(%1.2f,%1.2f), (sll,srr)=(%1.3f,%1.3f)\n",vsl,vsr,sll,srr);

            fl[0] = ql[dim+1];
            fl[1] = ql[1] * unl;
            fl[2] = ql[2] * unl;
            fl[3] = ql[3] * unl;
            fl[4] = (ql[4] + Pl) * unl;
            fl[dim + 1] += Pl;

            fr[0] = qr[dim+1];
            fr[1] = qr[1] * unr;
            fr[2] = qr[2] * unr;
            fr[3] = qr[3] * unr;
            fr[4] = (qr[4] + Pr) * unr;
            fr[dim + 1] += Pr;

            const double fl_mult = double(sll>0);
            const double fr_mult = double(srr<0);
            const double fhll_mult = double(srr>=0 and sll<=0);

            for (int j = 0; j < 5; ++j) {
                F[j] = (srr * fl[j] - sll * fr[j] + sll * srr * (qr[j] - ql[j])) / (srr - sll);

                F[j] = fhll_mult*F[j] + fl_mult*fl[j] + fr_mult*fr[j];
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
//        euler_solver::solve(_gamma, grid, q[idx], rhs[idx]);
        //euler_solver::solve_HLL(_gamma, variable_manager.dt(), grid, q[idx], rhs[idx]);
        euler_solver::solve_HLL_fast(_gamma, variable_manager.dt(), grid, q[idx], rhs[idx]);

    }

}

}
}
