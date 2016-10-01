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
datachunk::operator=(const datachunk & dataset)
{
    resize(dataset.range());
    pull(dataset);
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

}
}

