#include "rectilinear.h"
#include <fstream>
#include <cmath>

// QK includes
#include "lib/qkfunctions.h"

#define DIM_CHECK
#define INDEX_CHECK

namespace qk
{

namespace grid
{

rectilinear::rectilinear():
    qk::grid::grid()
{

}

rectilinear::rectilinear(const range & range, const double *startxs, const double *widths):
		qk::grid::grid(range,startxs,widths)
{
    _dxs.resize(num_dims(),0);
    for(int i = 0; i < _num_dims; i++){
        _dxs[i] = widths[i] / double(length(i));
    }
}

rectilinear::~rectilinear()
{

}

double
rectilinear::centroid(const int dim, const int index) const
{
#ifdef DIM_CHECK
    if(dim >= num_dims() || dim < 0){
        std::cout << "qk::grid::rectilinear::centroid : Requested dimension outside available\n";
        exit(EXIT_FAILURE);
    }
#endif
#ifdef INDEX_CHECK
    if(index >= upper(dim) || index < lower(dim)){
        std::cout << "qk::grid::rectilinear::centroid : Requested index outside available\n";
        exit(EXIT_FAILURE);
    }
#endif
    return _startxs[dim] + (index+0.5)*_dxs[dim];
}

double
rectilinear::dx(const int dim) const
{
#ifdef DIM_CHECK
    if(dim >= num_dims() || dim < 0){
        std::cout << "qk::grid::rectilinear::dx : Requested dimension outside available\n";
        exit(EXIT_FAILURE);
    }
#endif
    return _dxs[dim];
}

}
}

