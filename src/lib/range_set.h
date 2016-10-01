#ifndef _qk_lib_range_set_H
#define _qk_lib_range_set_H

// QK includes
#include "lib/indexer_interface.h"
#include "lib/range.h"

namespace qk
{

class range_set:
        public qk::indexer_interface<qk::range>
{
public:

    range_set();
    range_set(const qk::range & superRange, const int * numSubdivisions);

    ~range_set();

    void setup(const qk::range & superRange, const int * numSubdivisions);

private:
    qk::range _superRange;
};

}

#endif // _qk_lib_range_set_H
