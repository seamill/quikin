#include "face_basis.h"

// QK includes
#include "lib/exception.h"

namespace qk
{
namespace basis
{

face_basis::face_basis() :
        qk::basis::basis("face_basis", 0),
        _num_dims(0),
        _num_points_per_dim(0)
{

}

face_basis::~face_basis()
{
    // Nothing to do - data array is cleared by indexer interface
}

void face_basis::setup(const int num_dims, const int num_points_per_dim)
{
    _num_dims = num_dims;
    _num_points_per_dim = num_points_per_dim;
    _num_points = _num_dims * _num_points_per_dim;

    if (num_points_per_dim != 1) {
        throw qk::exception("qk::basis::face_basis::setup : Not yet setup for multiple points per dimensions.");
    }

    _normalized_face_points.push_back(0.);
    _normalized_quadrature_weights.push_back(1.);

}

}
}

