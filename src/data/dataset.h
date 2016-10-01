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

namespace qk
{
namespace data
{

class dataset:
        public qk::indexer_interface<extended_datachunk>
{
public:

    dataset();
    dataset(const qk::range & global_data_range);
    ~dataset();

    void resize(const qk::range & global_data_range);

    void sync();

//    void applyFunction( double (*func)(const QKIndexer & index) );

    //void writeVTK(const std::string & filename) const;

    //void clear();

protected:

    qk::range _data_range_global;
    qk::range _data_range_local;
    qk::range _chunk_range_local;

    std::vector<shared_datachunk> _lower_chunks;
    std::vector<shared_datachunk> _upper_chunks;

};

}
}

#endif // _qk_data_dataset_H
