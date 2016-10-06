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

    void resize(const qk::range & mesh_range, const qk::range & data_range, const int num_ghost_layers=2);

    void
    write_VTK(std::ofstream & file, const qk::range & range) const;

    int num_ghost_layers() const;

    const qk::range & internal_range() const
    {
        return _internal_range;
    }
    const qk::range & lower_internal_range(const int dim) const
    {
        return _lower_internal_ranges[dim];
    }
    const qk::range & upper_internal_range(const int dim) const
    {
        return _upper_internal_ranges[dim];
    }
    const qk::range & lower_external_range(const int dim) const
    {
        return _lower_external_ranges[dim];
    }
    const qk::range & upper_external_range(const int dim) const
    {
        return _upper_external_ranges[dim];
    }

protected:

    int _num_ghost_layers;

    qk::range _internal_range;

    std::vector<qk::range> _lower_internal_ranges;
    std::vector<qk::range> _lower_external_ranges;
    std::vector<qk::range> _upper_internal_ranges;
    std::vector<qk::range> _upper_external_ranges;

};

}
}

#endif // _qk_data_extended_datachunk_H
