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

    // Write data
    file << "LOOKUP_TABLE default" << std::endl;

    const double * d = _data.data();
    const int offset = range.lower(range.num_dims()-1);

    if(_mesh_range.num_dims()==3){
        int index[3];
        for(index[2]=range.lower(2); index[2] < range.upper(2); index[2]++){
            const int z_index = offset+(index[2] - lower(2)) * stride(2);
            for(index[1]=range.lower(1); index[1] < range.upper(1); index[1]++){
                const int yz_index = z_index + (index[1] - lower(1))*stride(1);
                for(index[0]=range.lower(0); index[0] < range.upper(0); index[0]++){
                    file << d[yz_index + (index[0] - lower(0))*stride(0)] << " ";
                }
            }
        }
    } else if(_mesh_range.num_dims()==2){
        int index[2];
        const int z_index = offset;
        for(index[1]=range.lower(1); index[1] < range.upper(1); index[1]++){
            const int yz_index = z_index + (index[1] - lower(1))*stride(1);
            for(index[0]=range.lower(0); index[0] < range.upper(0); index[0]++){
                file << d[yz_index + (index[0] - lower(0))*stride(0)] << " ";
            }
        }
    } else if(_mesh_range.num_dims()==1){
        int index[1];
        const int yz_index = offset;
        for(index[0]=range.lower(0); index[0] < range.upper(0); index[0]++){
            file << d[yz_index + (index[0] - lower(0))*stride(0)] << " ";
        }
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

void
datachunk::multiply_by_scalar(const qk::range & rng, const double scalar)
{
    for(qk::indexer idx = this->indexer(rng); idx.exists(); ++idx){
        (*this)[idx] *= scalar;
    }
}

}
}

