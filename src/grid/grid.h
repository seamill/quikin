#ifndef _qk_grid_grid_H
#define _qk_grid_grid_H

// STL includes
#include <fstream>

// QK includes
#include "lib/indexer.h"
#include "lib/range.h"
#include "basis/basis.h"

namespace qk
{
namespace grid
{

class grid: public qk::range
{
public:

    grid();

    grid(const range & range, const double *startxs, const double *widths);

    virtual ~grid();

    virtual void xc(const qk::indexer & idx, double * x) const = 0;

    virtual double start(const int dim) const;
    virtual double width(const int dim) const;

    virtual void write_vtk(std::ofstream & file, const qk::basis::basis & basis) const = 0;

    const qk::range & range() const
    {
        return *this;
    }

protected:

    void setup(const qk::range & range, const double *startxs, const double *widths);

    std::vector<double> _startxs;
    std::vector<double> _widths;

private:

};

}
}

#endif // _qk_grid_grid_H
