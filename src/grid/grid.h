#ifndef _qk_grid_grid_H
#define _qk_grid_grid_H

// QK includes
#include "lib/qkrange.h"

namespace qk
{
namespace grid
{

class grid: public range
{
public:

    grid();

    grid(const range & range, const double *startxs, const double *widths);

    virtual ~grid();

    virtual double start(const int dim) const;
    virtual double width(const int dim) const;

protected:

    void setup(const range & range, const double *startxs,
            const double *widths);

    std::vector<double> _startxs;
    std::vector<double> _widths;

private:

};

}
}

#endif // _qk_grid_grid_H
