#ifndef _qk_grid_structured_H
#define _qk_grid_structured_H

// QK includes
#include "lib/range.h"
#include "grid/grid.h"

namespace qk
{
namespace grid
{

class structured: public grid
{
public:

    structured();
    structured(const qk::range & range, const double *startxs, const double *widths);

    ~structured();

    void apply_bias_symmetric(const int dim, const double pow);
    void apply_bias_linear(const int dim, const double pow);

    void xc(const qk::indexer & idx, double * x) const
    {

    }

    double centroid(const int dim, const int index) const;
    double dx(const int dim, const int index) const;

protected:

    void bias_array(double * dxs, const int size, const int pow);
    void recalculate_xcs(const int dim);

    std::vector<std::vector<double> > _dxs;
    std::vector<std::vector<double> > _xcs;

private:

};

}
}

#endif // _qk_grid_structured_H
