#include <data/dataset.h>
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
    resize(qk::range());
}

dataset::~dataset()
{
    // Nothing to do here
}

void
dataset::sync()
{
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
dataset::resize(const qk::range & global_range)
{

    _data_range_global = global_range;

    if(global_range.num_dims() > 0){

        // For now this is a hack for a single process

        // Designate the local range belonging to this mpi rank
        _data_range_local = global_range;

        // Construct an mpi chunk size - only one chunk per process
        _chunk_range_local = global_range;
        for(int i = 0; i < _chunk_range_local.num_dims(); i++){
            _chunk_range_local.set(i,0,1);
        }

        qk::indexer_interface<extended_datachunk>::resize(_chunk_range_local);

        // Create local datachunk with numGhostLayers

        _data[0].resize(_data_range_local);

        // Create upper and lower shared datachunks (for mpi syncing)
    //    for(int i = 0; i < globalRange.getNumDims(); i++){

    //        QKRange lowerRange(_dataRange_local);
    //        lowerRange.set(i,lowerRange.getLower(i)-NUM_GHOST_LAYERS, lowerRange.getLower(i));
    //        _lowerChunks.push_back(QKSharedDatachunk(lowerRange));

    //        QKRange upperRange(_dataRange_local);
    //        upperRange.set(i, upperRange.getUpper(i), upperRange.getUpper(i)+NUM_GHOST_LAYERS);
    //        _upperChunks.push_back(QKSharedDatachunk(upperRange));

    //    }
    }
}


//void
//dataset::applyFunction( double (*func)(const QKIndexer & index) )
//{
//    for(QKIndexer indexer = getIndexer(); indexer.exists(); indexer.next()){
//        (*this)[indexer].applyFunction(func);
//    }
//}


}
}

