#include "qkrangeset.h"

// STL includes
#include <cmath>

namespace qk
{

range_set::range_set()
{
    // Do nothing
}

range_set::range_set(const qk::range & superRange, const int * numSubdivisions)
{
    setup(superRange, numSubdivisions);
}

range_set::~range_set()
{

}

void
range_set::setup(const qk::range &superRange, const int *numSubdivisions)
{
    // The goal is to break super range into a set of subdivisions in each dimension
    _superRange = superRange;
    resize(qk::range(_superRange.num_dims(), numSubdivisions));

    int subRangeSizePerDivision[_num_dims];
    for(int i = 0; i < _num_dims; i++){
        subRangeSizePerDivision[i] = int(std::ceil(float(superRange.length(i)) / float(numSubdivisions[i])));
    }

    int lowers[_num_dims];
    int uppers[_num_dims];

    for(qk::indexer indexer = this->indexer(); indexer.exists(); indexer.next()){
        for(int i = 0; i < _num_dims; i++){
            lowers[i] = _superRange.lower(i) + indexer[i] * subRangeSizePerDivision[i];
            uppers[i] = std::min(lowers[i] + subRangeSizePerDivision[i], _superRange.upper(i));
        }

        (*this)[indexer].resize(_num_dims, lowers, uppers);
    }

}

}
