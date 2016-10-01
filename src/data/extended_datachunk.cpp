#include "extended_datachunk.h"

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

extended_datachunk::extended_datachunk():
    datachunk()
{

}

extended_datachunk::~extended_datachunk()
{
    // Nothing to do - data array is cleared by indexer interface
}

void
extended_datachunk::resize(const qk::range & range)
{
    // While the interior range is what the dataset owns, the QKIndexInterface must contain ghost layers
    _internal_range = range;
    _lower_ranges.clear();
    _upper_ranges.clear();
    qk::range expanded_range = range;
    for(int i = 0; i < expanded_range.num_dims(); i++){
        expanded_range.expand(i, -NUM_GHOST_LAYERS, NUM_GHOST_LAYERS);
    }

    for(int i = 0; i < expanded_range.num_dims(); i++){
        qk::range lower_range(range);
        lower_range.set(i,lower_range.lower(i) - NUM_GHOST_LAYERS, lower_range.lower(i));
        _lower_ranges.push_back(lower_range);

        qk::range upper_range(range);
        upper_range.set(i, upper_range.upper(i), upper_range.upper(i) + NUM_GHOST_LAYERS);
        _upper_ranges.push_back(upper_range);
    }
    datachunk::resize(expanded_range);
}

int
extended_datachunk::num_ghost_layers()
{
    return NUM_GHOST_LAYERS;
}

}
}

