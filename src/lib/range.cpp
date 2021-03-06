#include "lib/range.h"

// STL includes
#include <cmath>
#include <algorithm>

// QK includes
#include "lib/functions.h"

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
    setup(range._num_dims, range._lower.data(), range._upper.data());
}

void
range::resize(const int num_dims, const int *lower, const int *upper)
{
    setup(num_dims, lower, upper);
}

bool
range::includes(const qk::range & rng) const
{
    if(rng.num_dims() != _num_dims){
        return false;
    }

    for(int i = 0;i<_num_dims;++i){
        const bool inclusive = (_lower[i] <= rng._lower[i]) and (_upper[i] >= rng._upper[i]);
        if(!inclusive){
            return false;
        }
    }
    return true;

}

void
range::setup(const int num_dims, const int *lower, const int *upper)
{
    // initialize
    _num_dims = num_dims;
    _lower.clear();
    _upper.clear();
    _stride.clear();

    // Setup new range

    if(_num_dims>0){
        _lower.resize(_num_dims,0);
        _upper.resize(_num_dims,0);
        _stride.resize(_num_dims,0);
        for(int i = 0; i < _num_dims; ++i){
            _lower[i] = lower[i];
            _upper[i] = upper[i];
        }
        _stride[_num_dims-1] = 1;
        for(int i = _num_dims-2; i >= 0; --i){
            _stride[i] = _stride[i+1] * length(i+1);
        }
    }

}

void
range::extrude(const int lower_extrude, const int upper_extrude)
{
    int num_dims = this->num_dims();
    std::vector<int> lower(num_dims+1,0), upper(num_dims+1,0);
    for(int i=0;i<num_dims;i++){
        lower[i] = this->lower(i);
        upper[i] = this->upper(i);
    }
    lower[num_dims] = lower_extrude;
    upper[num_dims] = upper_extrude;
    this->setup(num_dims+1,lower.data(), upper.data());
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
    if(_num_dims==0){return 0;}
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
range::operator=(const range & rng)
{
    setup(rng._num_dims, rng._lower.data(), rng._upper.data());
    return *this;
}

}

// rangeSet

