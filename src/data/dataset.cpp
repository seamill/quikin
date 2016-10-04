#include "data/dataset.h"

// STL includes
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>

#define NUM_GHOST_LAYERS 2

namespace qk
{
namespace data
{

dataset::dataset():
    qk::indexer_interface<extended_datachunk>()
{
    resize(qk::range(),qk::basis::basis(),1);
}

dataset::dataset(const qk::range & global_data_range, const qk::basis::basis & basis, const int num_components)
{
    resize(global_data_range, basis, num_components);
}

dataset::~dataset()
{
    // Nothing to do here
}

void
dataset::sync()
{
    throw qk::exception("qk::data::dataset::sync : Function really doesn't work yet.");

    // We iterate through all our datasets:
    // 1) copy ghost data from datasets
    // 2) send ghost data datasets
    // 3) recv ghost data datasets
    // 4) copy ghost data to datasets

    // We require:
    // 4 datasets per dimension per dataset (upper-send/upper-recv/lower-send/lower-recv)
    // 4 ranges per dimension per dataset (upper-send/upper-recv/lower-send/lower-recv)
}

void
dataset::resize(const qk::range & global_range, const qk::basis::basis & basis, const int num_components)
{
    _basis = basis;
    _mesh_range_global = global_range;
    _data_range = qk::range();
    _data_range.extrude(0,basis.num_points());
    _data_range.extrude(0,num_components);

    if(_mesh_range_global.volume() > 0){

        // For now this is a hack for a single process

        // Designate the local range belonging to this mpi rank
        _mesh_range_local = _mesh_range_global;

        // Construct an mpi chunk size - only one chunk per process
        _chunk_range_global = qk::range();
        for(int i = 0; i < _mesh_range_global.num_dims(); i++){
            _chunk_range_global.extrude(0,1);
        }

        qk::range chunk_range_local = _chunk_range_global;

        qk::indexer_interface<extended_datachunk>::resize(chunk_range_local);

        _data[0].resize(_mesh_range_local, _data_range);

    }
}

}
}

