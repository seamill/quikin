#ifndef QKRANGESET_H
#define QKRANGESET_H

// QK includes
#include "lib/qkrange.h"
#include "lib/qkindexerinterface.h"

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

#endif // QKRANGESET_H
