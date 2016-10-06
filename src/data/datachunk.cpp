#include "datachunk.h"

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

datachunk::datachunk():
    qk::indexer_interface<double>()
{

}

datachunk::~datachunk()
{
    // Nothing to do - data array is cleared by indexer interface
}

datachunk &
datachunk::operator=(const datachunk & chunk)
{
    resize(chunk._mesh_range, chunk._data_range);
    pull(chunk);
    return *this;
}

void
datachunk::write_VTK(std::ofstream & file, const qk::range & range) const
{

    file << "CELL_DATA " << range.volume() << std::endl;

    // Write data
    file << "SCALARS data float" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;

    // TODO: fix this so that it always does x,y,z
    for(qk::indexer idx = indexer(range); idx.exists(); idx.next()){
        file << (*this)[idx] << ' ';
    }
    file << '\n';

}

void
datachunk::resize(const qk::range & mesh_range, const qk::range & data_range)
{
    _mesh_range = mesh_range;
    _data_range = data_range;

    qk::range range = _mesh_range;
    for(int i = 0; i < _data_range.num_dims(); ++i){
        range.extrude(data_range.lower(i), data_range.upper(i));
    }

    qk::indexer_interface<double>::resize(range);

}

}
}

