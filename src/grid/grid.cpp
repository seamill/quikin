#include "grid.h"

// STL includes
#include <fstream>
#include <cmath>

// QK includes
#include "lib/functions.h"
#include "lib/exception.h"

namespace qk
{

namespace grid
{

grid::grid():
    qk::range()
{

}

grid::grid(const qk::range & range, const double *startxs, const double *widths):
    qk::range()
{
    setup(range, startxs, widths);
}

grid::~grid()
{

}

void
grid::setup(const qk::range & range, const double *startxs, const double *widths)
{
    resize(range);
    _startxs.resize(range.num_dims(),0);
    _widths.resize(range.num_dims(),0);
    for(int i = 0; i < _num_dims; i++){
        _widths[i] = widths[i];
        _startxs[i] = startxs[i];
        //_dxs[i] = widths[i] / double(getLength(i));
    }
}

double
grid::width(const int dim) const
{
#ifdef _QK_RANGE_CHECK_
    if(dim >= num_dims() || dim < 0){
        throw qk::exception("qk::grid::grid::width : Requested dimension out of range.");
    }
#endif
    return _widths[dim];
}


double
grid::start(const int dim) const
{
#ifdef _QK_RANGE_CHECK_
    if(dim >= num_dims() || dim < 0){
        throw qk::exception("qk::grid::grid::start : Requested dimension out of range.");
    }
#endif
    return _startxs[dim];
}

}
}

