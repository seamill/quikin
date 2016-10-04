#ifndef _qk_data_dataset_H
#define _qk_data_dataset_H

// STL includes
#include <string>

// QK includes
#include "lib/indexer.h"
#include "lib/indexer_interface.h"
#include "lib/range.h"
#include "data/extended_datachunk.h"
#include "data/shared_datachunk.h"
#include "basis/basis.h"

namespace qk
{
namespace data
{

class dataset:
        public qk::indexer_interface<extended_datachunk>
{
public:

    dataset();
    dataset(const qk::range & global_data_range, const qk::basis::basis & basis, const int num_components);
    ~dataset();

    void sync();

    const qk::basis::basis & basis() const {return _basis;}

protected:

    void resize(const qk::range & global_data_range, const qk::basis::basis & basis, const int num_components);

    qk::basis::basis _basis;
    qk::range _data_range;
    qk::range _mesh_range_global;
    qk::range _mesh_range_local;
    qk::range _chunk_range_global;

private:
    using qk::indexer_interface<extended_datachunk>::resize;

};

}
}

#endif // _qk_data_dataset_H
