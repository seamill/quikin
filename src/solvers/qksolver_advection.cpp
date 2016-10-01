#include "qksolver_advection.h"

// STL includes
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <chrono>
#include <iostream>
#include <stdlib.h>

// Vectorization includes
#include <immintrin.h>

// QK includes
#include "grid/rectilinear.h"
#include "lib/qkfunctions.h"

namespace qk
{
namespace solver
{

namespace advection_solver{

    void advection_iterate(const int start, const int end, const int offset, const double coeff, const double * __restrict__ n, double * __restrict__ np)
    {

        // Let compiler know of alignment
        // double is 8 bytes - AVX vectorizes 4 doubles -> 32 byte alignment required
        const double * __restrict__ a = (const double *) __builtin_assume_aligned((const void *) n, 32);
        double * __restrict__ b = (double *) __builtin_assume_aligned((void *) np, 32);

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


advection::advection()
{

}

advection::~advection()
{

}

void
advection::setup(const qk::grid::rectilinear & grid, const double *velocity)
{
    _grid = grid;

    const int numDims = grid.num_dims();
    _velocity.resize(numDims,0.0);
    for(int i = 0; i < numDims; i++){
        _velocity[i] = velocity[i];
    }

    _n.resize(grid);
    _ns.resize(grid);
    _np.resize(grid);

    // Now we initialize a gaussian in n

    const double dx[3] = {grid.dx(0), grid.dx(1), grid.dx(2)};
    const double xmin[3] = {grid.start(0), grid.start(1), grid.start(2)};
    const double xc[3] = {0., 0., 0.};
    const double sqrt2stdv = 0.2;

    double xs[3];
    double exponent;


    for(qk::indexer chunkIndexer = _n.indexer(); chunkIndexer.exists(); chunkIndexer.next()){
        qk::data::extended_datachunk & data = _n[chunkIndexer];
        for(qk::indexer indexer = data.indexer(); indexer.exists(); indexer.next()){
            xs[0] = (xmin[0] + dx[0] * indexer[0] - xc[0]) / sqrt2stdv;
            xs[1] = (xmin[1] + dx[1] * indexer[1] - xc[1]) / sqrt2stdv;
            xs[2] = (xmin[2] + dx[2] * indexer[2] - xc[2]) / sqrt2stdv;

            exponent = xs[0]*xs[0]+xs[1]*xs[1]+xs[2]*xs[2];

            data[indexer] = std::exp(-exponent);
        }
    }
}

void
advection::advance(const double time, const double dt)
{

    // Iterate through the domains in n
    for(qk::indexer indexer = _np.indexer(); indexer.exists(); indexer.next()){

        advection_solver::advection(
                    _grid,
                    _n[indexer].range(),
                    dt,
                    _velocity,
                    _n[indexer].data(),
                    _ns[indexer].data(),
                    _np[indexer].data()
                    );

        _n[indexer].swap(_np[indexer]);

    }
}

void
advection::write_VTK(const std::string & prefix, const std::string & suffix) const
{

    // Iterate through the domains in n
    for(qk::indexer indexer = _n.indexer(); indexer.exists(); indexer.next()){
        const qk::data::extended_datachunk & chunk = _n[indexer];
        const qk::range & range = chunk.internal_range();

        // Generate Filename

        // TODO: Add range information
        std::string rangeInfo = "";

        std::string filename = prefix + rangeInfo + suffix;

        std::ofstream file(filename.c_str());
        if(!file){
            std::cout << "QKDatachunk::writeCSV : File " << filename << " could not be created\n";
            exit(EXIT_FAILURE);
        }

        file << "# vtk DataFile Version 2.0\n";
        file << "Stuff for dataset\n";
        file << "ASCII\n";

        //_grid.writeVTK(file, range);
        chunk.write_VTK(file, range);

    }

}

}
}
