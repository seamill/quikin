#ifndef QKSOLVER_ADVECTION_H
#define QKSOLVER_ADVECTION_H

// QK includes
#include <data/dataset.h>
#include <solvers/solver.h>
#include "grid/rectilinear.h"

namespace qk
{

namespace solver
{

class advection:
        public solver
{
public:
	advection();
    ~advection();

    void setup(const qk::grid::rectilinear & grid, const double *velocity);

    void advance(const double time, const double dt);

    void write_VTK(const std::string & prefix, const std::string & suffix) const;

protected:

    //void updateDataset(const double time, const double dt, const QKDataset & n, QKDataset & ns, QKDataset & np) const;

    qk::data::dataset _n;
    qk::data::dataset _ns;
    qk::data::dataset _np;
    qk::grid::rectilinear _grid;
    std::vector<double> _velocity;
};

}
}

#endif // QKSOLVER_ADVECTION_H
