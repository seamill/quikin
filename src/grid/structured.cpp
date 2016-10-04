#include "structured.h"

// STL includes
#include <fstream>
#include <cmath>

// QK includes
#include "lib/functions.h"
#include "lib/exception.h"

#define DIM_CHECK
#define INDEX_CHECK

namespace qk
{

namespace grid
{

structured::structured():
		qk::grid::grid()
{

}

structured::structured(const qk::range & range, const double *startxs, const double *widths):
		qk::grid::grid(range,startxs,widths)
{
    _dxs.resize(num_dims());
    _xcs.resize(num_dims());
    for(int i = 0; i < _num_dims; i++){
        const double dx = widths[i] / double(length(i));
        _dxs[i].resize(range.length(i),dx);
        _xcs[i].resize(range.length(i),0.);
        recalculate_xcs(i);
    }

}

structured::~structured()
{

}

void
structured::bias_array(double * dxs, const int size, const int pow)
{
    // Thi
    double xs[size+1];
    xs[0] = 0.;
    const double dV = 1.0 / double(size);
    const double invPow = 1. / double(pow);
    for(int i = 0; i < size; i++){
        xs[i+1]=std::pow(std::pow(xs[i],pow) + dV, invPow);
        dxs[i] = xs[i+1]-xs[i];
    }
}

void
structured::apply_bias_symmetric(const int dim, const double pow)
{
#ifdef DIM_CHECK
    if(dim >= num_dims() || dim < 0){
        throw qk::exception("qk::grid::structured::apply_bias_symmetric : Requested dimension out of range.");
    }
#endif
    const int L = width(dim)/2.;
    const int N = length(dim);

    if(N %2 == 1){
        throw qk::exception("qk::grid::structured::apply_bias_symmetric : Can only bias even grids.");
    }

    const int n = N/2;
    double n_dxs[n];
    bias_array(n_dxs,n,pow);
    double * dxs = &(_dxs[dim][0]);
    for(int i = 0; i < n; i++){
        dxs[i] = L*n_dxs[i];
        dxs[N-1-i] = dxs[i];
    }

    recalculate_xcs(dim);
}

void
structured::apply_bias_linear(const int dim, const double pow)
{
#ifdef DIM_CHECK
    if(dim >= num_dims() || dim < 0){
        throw qk::exception("qk::grid::structured::apply_bias_linear : Requested dimension out of range.");
    }
#endif
    const int L = width(dim);
    const int N = length(dim);
    double n_dxs[N];
    bias_array(n_dxs,N,pow);
    double * dxs = &(_dxs[dim][0]);
    for(int i = 0; i < N; i++){
        dxs[i] = L*n_dxs[N-1-i];
    }

    recalculate_xcs(dim);

}

double
structured::centroid(const int dim, const int index) const
{
#ifdef DIM_CHECK
    if(dim >= num_dims() || dim < 0){
        throw qk::exception("qk::grid::structured::centroid : Requested dimension out of range.");
    }
#endif
#ifdef INDEX_CHECK
    if(index >= upper(dim) || index < lower(dim)){
        throw qk::exception("qk::grid::structured::centroid : Requested index out of range.");
    }
#endif
    return _xcs[dim][index];
}

double
structured::dx(const int dim, const int index) const
{
#ifdef DIM_CHECK
    if(dim >= num_dims() || dim < 0){
        throw qk::exception("qk::grid::structured::dx : Requested dimension out of range.");
    }
#endif
#ifdef INDEX_CHECK
    if(index >= upper(dim) || index < lower(dim)){
        throw qk::exception("qk::grid::structured::dx : Requested index out of range.");
    }
#endif
    return _dxs[dim][index];
}

void
structured::recalculate_xcs(const int dim)
{
#ifdef DIM_CHECK
    if(dim >= num_dims() || dim < 0){
        throw qk::exception("qk::grid::structured::recalculate_xcs : Requested dimension out of range.");
    }
#endif
    const int N = length(dim);
    const double * dxs = &(_dxs[dim][0]);
    double * xcs = &(_xcs[dim][0]);

    double pxc = _startxs[dim] + 0.5*dxs[0];
    xcs[0] = pxc;
    for(int i = 1; i < N; i++){
        pxc += 0.5*(dxs[i-1] + dxs[i]);
        xcs[i] = pxc;
    }
}

}
}

