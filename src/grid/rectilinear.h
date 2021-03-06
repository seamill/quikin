#ifndef _qk_grid_rectilinear_H
#define _qk_grid_rectilinear_H

// QK includes
#include "lib/range.h"
#include "grid/grid.h"

namespace qk
{
namespace grid
{

class rectilinear: public grid
{
public:

    rectilinear();
    rectilinear(const qk::range & range, const double *startxs, const double *widths);

    ~rectilinear();

    void xc(const qk::indexer & idx, double * x) const;

    double centroid(const int dim, const int index) const;

    double dx(const int dim) const;

    void write_vtk(std::ofstream & file, const qk::basis::basis & basis) const;

protected:

    std::vector<double> _dxs;

private:

};

}
}

#endif // _qk_grid_rectilinear_H
