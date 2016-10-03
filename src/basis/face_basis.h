#ifndef _qk_basis_face_basis_H
#define _qk_basis_face_basis_H

// STL includes
#include <vector>

// QK includes
#include "basis/basis.h"

namespace qk
{
namespace basis
{

class face_basis: public qk::basis::basis
{
public:

    face_basis();
    virtual ~face_basis();

    void setup(const int num_dims, const int num_points_per_dim);

protected:

    int _num_dims;
    int _num_points_per_dim;
    std::vector<double> _normalized_face_points;
    std::vector<double> _normalized_quadrature_weights;

};

}
}

#endif // _qk_basis_face_basis_H
