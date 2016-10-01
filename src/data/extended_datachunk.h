#ifndef _qk_data_extended_datachunk_H
#define _qk_data_extended_datachunk_H

// STL includes
#include <string>

// QK includes
#include "lib/indexer.h"
#include "lib/indexer_interface.h"
#include "lib/range.h"
#include "data/datachunk.h"

namespace qk
{
namespace data
{

class extended_datachunk: public datachunk
{
public:

    extended_datachunk();
    ~extended_datachunk();

    void resize(const qk::range & internal_range);

    static int num_ghost_layers();

    const qk::range & internal_range() const
    {
        return _internal_range;
    }
    const qk::range & lower_ghost_range(const int dim) const
    {
        return _lower_ranges[dim];
    }
    const qk::range & upper_ghost_range(const int dim) const
    {
        return _upper_ranges[dim];
    }

protected:

    qk::range _internal_range;

    std::vector<qk::range> _lower_ranges;
    std::vector<qk::range> _upper_ranges;

};

}
}

#endif // _qk_data_extended_datachunk_H
