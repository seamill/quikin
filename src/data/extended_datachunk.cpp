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
extended_datachunk::resize(const qk::range & mesh_range, const qk::range & data_range)
{


    qk::range extended_mesh_range = mesh_range;
    for(int i = 0; i < extended_mesh_range.num_dims(); i++){
        extended_mesh_range.expand(i, -NUM_GHOST_LAYERS, NUM_GHOST_LAYERS);
    }

    qk::data::datachunk::resize(extended_mesh_range, data_range);

    const qk::range & range = this->range();

    // While the interior range is what the dataset owns, the QKIndexInterface must contain ghost layers

    _internal_range = qk::range();
    for(int i = 0; i < mesh_range.num_dims(); i++){
        _internal_range.extrude(mesh_range.lower(i), mesh_range.upper(i));
    }
    for(int i = 0; i < data_range.num_dims(); i++){
        _internal_range.extrude(data_range.lower(i), data_range.upper(i));
    }

    _lower_internal_ranges.clear();
    _lower_external_ranges.clear();
    _upper_internal_ranges.clear();
    _upper_external_ranges.clear();

    for(int i = 0; i < extended_mesh_range.num_dims(); i++){
        qk::range lower_external_range(_internal_range);
        lower_external_range.set(i, lower_external_range.lower(i) - num_ghost_layers(), lower_external_range.lower(i));
        _lower_external_ranges.push_back(lower_external_range);

        qk::range lower_internal_range(_internal_range);
        lower_internal_range.set(i, lower_internal_range.lower(i), lower_internal_range.lower(i) + num_ghost_layers());
        _lower_internal_ranges.push_back(lower_internal_range);

        qk::range upper_internal_range(_internal_range);
        upper_internal_range.set(i, upper_internal_range.upper(i) - num_ghost_layers(), upper_internal_range.upper(i));
        _upper_internal_ranges.push_back(upper_internal_range);

        qk::range upper_external_range(_internal_range);
        upper_external_range.set(i, upper_external_range.upper(i), upper_external_range.upper(i) + num_ghost_layers());
        _upper_external_ranges.push_back(upper_external_range);
    }

}

int
extended_datachunk::num_ghost_layers()
{
    return NUM_GHOST_LAYERS;
}

}
}

