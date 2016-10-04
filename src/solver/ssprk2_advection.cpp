#include "ssprk2_advection.h"

// STL include
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <chrono>
#include <iostream>
#include <stdlib.h>

// QK includes
#include "lib/functions.h"
#include "lib/exception.h"
#include "grid/rectilinear.h"
#include "variable/variable_id.h"
#include "variable/variable.h"
#include "data/functions.h"

namespace qk
{
namespace solver
{

namespace advection_solver{

    void advection_iterate(const int start, const int end, const int offset, const double coeff, const double * __restrict__ n, double * __restrict__ np)
    {

        // Let compiler know of alignment
        // double is 8 bytes - AVX vectorizes 4 doubles -> 32 byte alignment required
        const double * __restrict__ a = n;//(const double *) __builtin_assume_aligned((const void *) n, 32);
        double * __restrict__ b = np;//(double *) __builtin_assume_aligned((void *) np, 32);

        // Iterate - should autovectorize
#pragma omp parallel for
        for(int i = start; i < end; i++){
            // Center cell value
            const double vc = a[i-offset];
            // Left cell value
            const double vn = a[i-2*offset];
            // Right cell value
            const double vp = a[i];

            const double r = (vp != vc) ? (vc - vn) / (vp - vc) : 1.0; // 2 add, 1 div : 3 FLO

            //const double phi = std::max(0.0, std::min(2. * r, std::min(2., (2. + r) / 3.))); // 1 add, 1 mul, 1 div, 3 min/max : 4 FLO
            const double phi = r*(r + 1.0) / (r*r + 1.0);
            const double dflux = coeff * (vc + 0.25 * phi * (vp - vn)); // 2 add, 3 mul : 5 FLO

//            if(phi != phi){
//                printf("(% e, % e, % e) : phi % e, Flux % e\n", vn, vc, vp, phi, dflux);
//            }

            b[i-offset] -= dflux; // 1 add : 1 FLO
            b[i] += dflux; // 1 add : 1FLO

        }
    }

    void advection(
            const qk::grid::rectilinear & grid,
            const qk::range & totalRange,
            const double & dt,
            const std::vector<double> & velocity,
            const double * __restrict__ n_data,
            double * __restrict__ ns_data,
            double * __restrict__ np_data
            )
    {
        const int numDims = totalRange.num_dims();
        const int volume = totalRange.volume();
        int stride[numDims];
        stride[0] = 1;
        for(int i = 1; i < numDims; i++){
            stride[i] = stride[i-1] * totalRange.length(i-1);
        }

        const double * n = (const double *) n_data;//__builtin_assume_aligned(n_data, 32);
        double * ns = (double *) ns_data;//__builtin_assume_aligned(ns_data, 32);
        double * np = (double *) np_data;//__builtin_assume_aligned(np_data, 32);

        // Initialize ns - maybe use memcpy?
        //memcpy(ns, n, volume * sizeof(double));
        for(int i = 0; i < volume; i++){
            ns[i] = n[i];
            //printf("ns = % e\n", ns[i]);
        }

        // Apply fluxes
        for(int dim = 0; dim < numDims; dim++){
            const int offset = ((velocity[dim] > 0) ? 1 : -1) * stride[dim];
            const int start = std::max(0, 2 * offset);
            const int end = std::min(volume, volume-2*offset);
            advection_iterate(start, end, offset, dt * velocity[dim] / grid.dx(dim), n, ns);
        }

        // Initialize np
        for(int i = 0; i < volume; i++){
            np[i] = 0.5*(n[i]+ns[i]);
            //printf("np = % e\n", np[i]);
        }

        // Apply fluxes
        for(int dim = 0; dim < numDims; dim++){
            const int offset = ((velocity[dim] > 0) ? 1 : -1) * stride[dim];
            const int start = std::max(0, 2 * offset);
            const int end = std::min(volume, volume-2*offset);
            advection_iterate(start, end, offset, 0.5 * dt * velocity[dim] / grid.dx(dim), ns, np);
        }

    }

}


ssprk2_advection::ssprk2_advection()
{

}

ssprk2_advection::~ssprk2_advection()
{

}

void
ssprk2_advection::setup(const std::vector<double> & velocity)
{
    _velocity = velocity;
}

void
ssprk2_advection::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{
    if(_input_variable_ids.size() != 1){
        throw qk::exception("qk::solver::ssprk2_advection::solve : Input must contain one variable.");
    }

    if(_output_variable_ids.size() != 2){
        throw qk::exception("qk::solver::ssprk2_advection::solve : Output must contain two variables.");
    }

    const qk::grid::rectilinear & grid = dynamic_cast<const qk::grid::rectilinear &>(variable_manager.grid());

    const double time = variable_manager.time();
    const double dt = variable_manager.dt();

    const qk::variable::variable_id & input_id = _input_variable_ids[0];
    const qk::variable::variable & n_0 = variable_manager.input_variable(input_id);

    const qk::variable::variable_id & partial_id = _output_variable_ids[0];
    qk::variable::variable & n_1 = variable_manager.output_variable(partial_id);

    const qk::variable::variable_id & output_id = _output_variable_ids[1];
    qk::variable::variable & n_2 = variable_manager.output_variable(output_id);

    if(tag == STAGE_0){
        for(qk::indexer chunk_idx = n_0.indexer(); chunk_idx.exists(); chunk_idx.next()){
            const qk::data::extended_datachunk & n_0_chunk = n_0[chunk_idx];
            qk::data::extended_datachunk & n_1_chunk = n_1[chunk_idx];

//            stage_0(dt, grid, n_0_chunk, n_1_chunk);

            qk::data::extended_datachunk & n_2_chunk = n_2[chunk_idx];
            stage_0(dt, grid, n_0_chunk, n_2_chunk);
        }
    } else if (tag == STAGE_1){

        return;

        for(qk::indexer chunk_idx = n_0.indexer(); chunk_idx.exists(); chunk_idx.next()){
            const qk::data::extended_datachunk & n_0_chunk = n_0[chunk_idx];
            const qk::data::extended_datachunk & n_1_chunk = n_1[chunk_idx];
            qk::data::extended_datachunk & n_2_chunk = n_2[chunk_idx];

            stage_1(dt, grid, n_0_chunk, n_1_chunk, n_2_chunk);
        }
    }

}

void
ssprk2_advection::stage_0(const double dt, const qk::grid::rectilinear & grid, const qk::data::extended_datachunk & n0, qk::data::extended_datachunk & n1) const
{
    const int num_dims = grid.num_dims();
    int stride[num_dims];
    for(int i = 0; i < num_dims; i++){
        stride[i] = n0.range().stride(i);
    }

    const double * n = (const double *) n0.data();//__builtin_assume_aligned(n_data, 32);
    double * ns = (double *) n1.data();//__builtin_assume_aligned(ns_data, 32);

    // Initialize n1 with n0
    qk::data::full_copy(n0,n1);

    // Apply fluxes
    for(int dim = 0; dim < num_dims; dim++){
        const int offset = ((_velocity[dim] > 0) ? 1 : -1) * stride[dim];
        const int start = std::max(0, n0.num_ghost_layers() * offset);
        const int end = std::min(n0.volume(), n0.volume()-n0.num_ghost_layers()*offset);
        advection_solver::advection_iterate(start, end, offset, dt * _velocity[dim] / grid.dx(dim), n, ns);
    }
}

void
ssprk2_advection::stage_1(const double dt, const qk::grid::rectilinear & grid, const qk::data::extended_datachunk & n0, const qk::data::extended_datachunk & n1, qk::data::extended_datachunk & n2) const
{

    const int num_dims = grid.num_dims();
    const int volume = n0.volume();
    int stride[num_dims];
    for(int i = 0; i < num_dims; i++){
        stride[i] = n0.range().stride(i);
    }

    const double * n = (const double *) n0.data();
    const double * ns = (const double *) n1.data();
    double * np = (double *) n2.data();

    // Initialize np
    for(int i = 0; i < volume; i++){
        np[i] = 0.5*(n[i]+ns[i]);
    }

    // Apply fluxes
    for(int dim = 0; dim < num_dims; dim++){
        const int offset = ((_velocity[dim] > 0) ? 1 : -1) * stride[dim];
        const int start = std::max(0, n0.num_ghost_layers() * offset);
        const int end = std::min(volume, volume-n0.num_ghost_layers()*offset);
        advection_solver::advection_iterate(start, end, offset, 0.5 * dt * _velocity[dim] / grid.dx(dim), ns, np);
    }

}

}
}
