#include "rectilinear.h"

// STL includes
#include <fstream>
#include <cmath>

// QK includes
#include "lib/functions.h"
#include "lib/exception.h"
#include "basis/volume_average.h"

namespace qk
{

namespace grid
{

rectilinear::rectilinear():
    qk::grid::grid()
{

}

rectilinear::rectilinear(const qk::range & range, const double *startxs, const double *widths):
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

void
rectilinear::write_vtk(std::ofstream & file, const qk::basis::basis & basis) const
{

    if(basis == qk::basis::volume_average()){

        int nx=1, ny=1, nz=1;
        if(_num_dims > 0){nx = length(0)+1;}
        if(_num_dims > 1){ny = length(1)+1;}
        if(_num_dims > 2){nz = length(2)+1;}

        double cx=0, cy=0, cz=0;
        if(_num_dims > 0){cx = _startxs[0] + _widths[0] / 2.;}
        if(_num_dims > 1){cy = _startxs[1] + _widths[1] / 2.;}
        if(_num_dims > 2){cz = _startxs[2] + _widths[2] / 2.;}

        double sx=0, sy=0, sz=0;
        if(_num_dims > 0){sx = _dxs[0];}
        if(_num_dims > 1){sy = _dxs[1];}else{sy=sx;}
        if(_num_dims > 2){sz = _dxs[2];}else{sz=sx;}

        file << "DATASET STRUCTURED_POINTS\n";
        file << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
        file << "ORIGIN " << cx << " " << cy << " " << cz << "\n";
        file << "SPACING " << sx << " " << sy << " " << sz << "\n";


    } else {
        throw qk::exception("qk::grid::rectilinear::write_vtk : Basis not allowed (i.e. not setup yet).");
    }

}

void
rectilinear::xc(const qk::indexer & idx, double * x) const
{
    for(int i=0;i<this->num_dims();i++){
        x[i] = _startxs[i] + (idx.index(i)+0.5)*_dxs[i];
    }
}

double
rectilinear::centroid(const int dim, const int index) const
{
#ifdef _QK_RANGE_CHECK_
    if(dim >= num_dims() || dim < 0){
        throw qk::exception("qk::grid::rectilinear::centroid : Requested dimension out of range.");
    }
    if(index >= upper(dim) || index < lower(dim)){
        throw qk::exception("qk::grid::rectilinear::centroid : Requested index out of range.");
    }
#endif
    return _startxs[dim] + (index+0.5)*_dxs[dim];
}

double
rectilinear::dx(const int dim) const
{
#ifdef _QK_RANGE_CHECK_
    if(dim >= num_dims() || dim < 0){
        throw qk::exception("qk::grid::rectilinear::dx : Requested dimension out of range.");
    }
#endif
    return _dxs[dim];
}

}
}

