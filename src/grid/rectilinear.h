#ifndef _qk_grid_rectilinear_H
#define _qk_grid_rectilinear_H

// QK includes
#include "grid/grid.h"
#include "lib/qkrange.h"

namespace qk
{
namespace grid
{

class rectilinear: public grid
{
public:

    rectilinear();
    rectilinear(const range & range, const double *startxs,
            const double *widths);

    ~rectilinear();

    double centroid(const int dim, const int index) const;

    double dx(const int dim) const;

protected:

    std::vector<double> _dxs;

private:

};

}
}

#endif // _qk_grid_rectilinear_H
