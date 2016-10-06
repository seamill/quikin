#include "extended_datachunk.h"

// STL includes
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>

namespace qk
{
namespace data
{

extended_datachunk::extended_datachunk():
    datachunk(),
    _num_ghost_layers(2)
{

}

extended_datachunk::~extended_datachunk()
{
    // Nothing to do - data array is cleared by indexer interface
}

void
extended_datachunk::resize(const qk::range & mesh_range, const qk::range & data_range, const int num_ghost_layers)
{
    _num_ghost_layers = num_ghost_layers;

    qk::range extended_mesh_range = mesh_range;
    for(int i = 0; i < extended_mesh_range.num_dims(); i++){
        extended_mesh_range.expand(i, -this->num_ghost_layers(), this->num_ghost_layers());
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
        lower_external_range.set(i, lower_external_range.lower(i) - this->num_ghost_layers(), lower_external_range.lower(i));
        _lower_external_ranges.push_back(lower_external_range);

        qk::range lower_internal_range(_internal_range);
        lower_internal_range.set(i, lower_internal_range.lower(i), lower_internal_range.lower(i) + this->num_ghost_layers());
        _lower_internal_ranges.push_back(lower_internal_range);

        qk::range upper_internal_range(_internal_range);
        upper_internal_range.set(i, upper_internal_range.upper(i) - this->num_ghost_layers(), upper_internal_range.upper(i));
        _upper_internal_ranges.push_back(upper_internal_range);

        qk::range upper_external_range(_internal_range);
        upper_external_range.set(i, upper_external_range.upper(i), upper_external_range.upper(i) + this->num_ghost_layers());
        _upper_external_ranges.push_back(upper_external_range);
    }

}

void
extended_datachunk::write_VTK(std::ofstream & file, const qk::range & range) const
{

    qk::data::datachunk::write_VTK(file,range);

//    file << "CELL_DATA " << _internal_range.volume() << std::endl;
//
//    // Write data
//    file << "SCALARS data float" << std::endl;
//    file << "LOOKUP_TABLE default" << std::endl;
//
//    // TODO: fix this so that it always does x,y,z
//    for(qk::indexer idx = indexer(_internal_range); idx.exists(); idx.next()){
//        file << (*this)[idx] << ' ';
//    }
//    file << '\n';

}

int
extended_datachunk::num_ghost_layers() const
{
    return _num_ghost_layers;
}

}
}

