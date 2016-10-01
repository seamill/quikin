#include "qkrange.h"

// STL includes
#include <cmath>
#include <algorithm>

// QK includes
#include "lib/qkfunctions.h"

// range

namespace qk
{

range::range():
    _num_dims(0)
{

}

range::range(const range & range):
    _num_dims(0)
{
    resize(range);
}

range::range(const int num_dims, const int * dims):
    _num_dims(0)
{
    int lower[num_dims];
    std::fill_n(lower,num_dims,0);
    setup(num_dims, lower, dims);
}

range::range(const int num_dims, const int * lower, const int * upper):
    _num_dims(0)
{
    setup(num_dims, lower, upper);
}

range::~range()
{

}

void
range::resize(const range &range)
{
    setup(range._num_dims, &(range._lower[0]), &(range._upper[0]));
}

void
range::resize(const int num_dims, const int *lower, const int *upper)
{
    setup(num_dims, lower, upper);
}

void
range::setup(const int num_dims, const int *lower, const int *upper)
{
    // Clear range
    _lower.clear();
    _upper.clear();
    _stride.clear();

    // Setup new range
    _num_dims = num_dims;
    if(_num_dims){
        _lower.resize(_num_dims);
        _upper.resize(_num_dims);
        _stride.resize(_num_dims);
        for(int i = 0; i < _num_dims; i++){
            _lower[i] = lower[i];
            _upper[i] = upper[i];
        }
        _stride[_num_dims-1] = 1;
        for(int i = _num_dims-2; i >= 0; i--){
            _stride[i] = _stride[i+1] * length(i+1);
        }
    }

}

void
range::expand(const int dim, const int lowerExpand, const int upperExpand)
{
    set(dim, _lower[dim]+lowerExpand, _upper[dim]+upperExpand);
}

void
range::shift(const int dim, const int shift)
{
    set(dim, _lower[dim]+shift, _upper[dim]+shift);
}

int
range::num_dims() const
{
    return _num_dims;
}

int
range::length(const int i) const
{
    return _upper[i] - _lower[i];
}

int
range::stride(const int i) const
{
    return _stride[i];
}

int
range::volume() const
{
    if(!_num_dims){return 0;}
    return _stride[0] * length(0);

}

void
range::set(const int dim, const int new_lower, const int new_upper)
{
    int lower[_num_dims], upper[_num_dims];
    for(int i = 0; i < _num_dims; i++){
        lower[i] = _lower[i];
        upper[i] = _upper[i];
    }
    lower[dim] = new_lower;
    upper[dim] = new_upper;

    setup(_num_dims, lower, upper);
}

range &
range::operator=(const range & range)
{
    resize(range);
    return *this;
}

}

// rangeSet

