#ifndef QKSOLVER_H
#define QKSOLVER_H

// QK includes
#include <data/dataset.h>
#include <grid/grid.h>

namespace qk
{
namespace solver
{

class solver
{
public:
    solver();

    virtual ~solver();

    virtual void advance(const double time, const double dt) = 0;

    virtual void write_VTK(const std::string & prefix, const std::string & suffix) const = 0;

protected:

};

}
}

#endif // QKSOLVER_H
