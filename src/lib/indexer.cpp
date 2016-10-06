#include "indexer.h"

// QK includes
#include "lib/exception.h"

namespace qk
{

indexer::indexer() :
        qk::range(),
        _sub_range(),
        _linear_index(0),
        _exists(false)
{

}

indexer::indexer(const qk::range & sub_range, const qk::range & sup_range) :
        qk::range(sup_range),
        _sub_range(sub_range),
        _linear_index(0),
        _exists(false)
{
    _index.resize(_num_dims,0);

    if (_num_dims > 0) {
        _exists = this->includes(_sub_range);
        if(_exists){
            for(int i = 0; i < _num_dims; i++){
                _index[i] = sub_range.lower(i);
            }
            _linear_index = 0;
            for (int i = 0; i < _num_dims; ++i){
                _linear_index += (_index[i] - _lower[i]) * _stride[i];
            }
        }
    }
}

indexer::~indexer()
{

}

const int &
indexer::operator[](const int & dim) const
{
#if _QK_RANGE_CHECK_
    if(dim < 0 or dim > _num_dims){
        throw qk::exception("qk::indexer::operator[] : Dimension out of range.");
    }
#endif
    return _index[dim];
}

int indexer::index(const int dim) const
{
#if _QK_RANGE_CHECK_
    if(dim < 0 or dim > _num_dims){
        throw qk::exception("qk::indexer::index : Dimension out of range.");
    }
#endif
    return _index[dim];
}

indexer &
indexer::operator=(const indexer & indexer)
{
    resize(dynamic_cast<const qk::range &>(indexer));
    _sub_range = indexer._sub_range;
    _exists = indexer._exists;
    _index.resize(_num_dims);
    _stride.resize(_num_dims);
    for (int i = _num_dims - 1; i >= 0; i--) {
        _index[i] = indexer._index[i];
        _stride[i] = indexer._stride[i];
    }
    _linear_index = indexer._linear_index;
    return *this;
}

indexer &
indexer::operator++()
{
    this->next();
    return *this;
}

void indexer::increment(const int dim)
{
#if _QK_RANGE_CHECK_
    if(dim < 0 or dim > _num_dims){
        throw qk::exception("qk::indexer::increment : Dimension out of range.");
    }
#endif
    _linear_index += _stride[dim];
    _index[dim]++;
}

void indexer::decrement(const int dim)
{
#if _QK_RANGE_CHECK_
    if(dim < 0 or dim > _num_dims){
        throw qk::exception("qk::indexer::decrement : Dimension out of range.");
    }
#endif
    _linear_index -= _stride[dim];
    _index[dim]--;
}

void indexer::next()
{
    // Indexer moves along outermost row and increments internal rows as needed
    recursive_increment_dim(_num_dims - 1);
}

void indexer::recursive_increment_dim(const int & dim)
{
    if (dim < 0) {
        _exists = false;
        return;
    }
    if (_index[dim] == _sub_range.upper(dim) - 1) {
        _index[dim] = _sub_range.lower(dim);
        _linear_index -= (_sub_range.length(dim) - 1) * _stride[dim];
        recursive_increment_dim(dim - 1);
        return;
    }
    increment(dim);
}

int indexer::linear_index() const
{
    return _linear_index;
}

bool indexer::exists() const
{
    return _exists;
}

}
