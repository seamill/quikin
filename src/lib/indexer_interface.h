#ifndef _qk_lib_indexer_interface_H
#define _qk_lib_indexer_interface_H

// QK includes
#include "lib/indexer.h"
#include "lib/range.h"
#include "lib/exception.h"

namespace qk
{

template<typename T>
class indexer_interface: public qk::range
{
public:

    indexer_interface(const qk::range & range = qk::range())
    {
        resize(range);
    }
    indexer_interface(const indexer_interface & interface)
    {
        resize(interface.range());
    }

    virtual ~indexer_interface()
    {

    }

    void swap(indexer_interface & interface)
    {
        //TODO: This check needs to be a bit more rugged.. doesn't check if internal elements are swappable
        if (this->range().volume() != interface.range().volume()) {
            throw qk::exception("qk::indexer_interface::swap : Volume mismatch.");
        }
        _data.swap(interface._data);
    }

    void pull(const indexer_interface<T> & interface,
            const qk::range & from_range, const qk::range & to_range)
    {
        if (from_range.volume() != to_range.volume()) {
            throw qk::exception("qk::indexer_interface::pull : Volume mismatch.");
        }

        // TODO: this needs to be sped up
        qk::indexer src_indexer = interface.indexer(from_range);
        qk::indexer dst_indexer = indexer(to_range);
        while (src_indexer.exists()) {
            (*this)[dst_indexer] = interface[src_indexer];
            dst_indexer.next();
            src_indexer.next();
        }
    }

    void pull(const indexer_interface<T> & interface)
    {
        // We assume interface range has the same volume as range
        // TODO: this needs to be sped up
        // Cant use memcopy if T has dynamic memory
        qk::indexer src_indexer = interface.indexer();
        qk::indexer dst_indexer = indexer();
        while (src_indexer.exists()) {
            (*this)[dst_indexer] = interface[src_indexer];
            dst_indexer.next();
            src_indexer.next();
        }
    }

    virtual void resize(const qk::range & range)
    {
        qk::range::resize(range);
        _data.resize(volume());
    }

    const qk::range & range() const
    {
        return dynamic_cast<const qk::range &>(*this);
    }

    qk::indexer indexer(const qk::range & range) const
    {
        return qk::indexer(range, this->range());
    }

    qk::indexer indexer() const
    {
        return qk::indexer(this->range(), this->range());
    }

    const T &
    operator [](const qk::indexer & index) const
    {
        return _data[index.linear_index()];
    }

    T &
    operator [](const qk::indexer & index)
    {
        return _data[index.linear_index()];
    }

    const T * data() const
    {
        return &(_data[0]);
    }
    T * data()
    {
        return &(_data[0]);
    }

    const T * data(const int * index) const
    {
        // NOTE : index is in global coordinates
        int idx = 0;
        for (int i = 0; i < _num_dims; i++) {
            idx += (index[i] - _lower[i]) * _stride[i];
        }
        return &(_data[0]) + idx;
    }

    T *
    data(const int * index)
    {
        // NOTE : index is in global coordinates
        int idx = 0;
        for (int i = 0; i < _num_dims; i++) {
            idx += (index[i] - _lower[i]) * _stride[i];
        }
        return &(_data[0]) + idx;
    }

protected:

    std::vector<T> _data;

};

}

#endif // _qk_lib_indexer_interface_H
