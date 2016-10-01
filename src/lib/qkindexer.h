#ifndef QKINDEXER_H
#define QKINDEXER_H

// QK incldues
#include "lib/qkrange.h"

namespace qk
{

class indexer:
        public qk::range
{
public:
    // Indexers are used to iterate around a subRange within a superRange
	indexer();
	indexer(const qk::range & sub_range, const qk::range & sup_range);

    ~indexer();

    void increment(const int dim);
    void decrement(const int dim);

    void next();
    bool exists() const;

    int linear_index() const;

    const int & operator[](const int & dim) const;
    indexer & operator=(const indexer & idx);

protected:

    void recursive_increment_dim(const int & dim);

    qk::range _sub_range;

    int _linear_index;
    bool _exists;
    std::vector<int> _index;
    std::vector<int> _stride;


};

}

#endif // QKINDEXER_H