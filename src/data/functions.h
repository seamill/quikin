#ifndef _qk_data_functions_H
#define _qk_data_functions_H

// STL includes
#include <string>

// QK includes
#include "lib/indexer.h"
#include "lib/indexer_interface.h"
#include "lib/range.h"

namespace qk
{
namespace data
{

    inline void
    copy_to(const qk::data::datachunk & from_chunk, const qk::range & from_range, qk::data::datachunk & to_chunk, const qk::range & to_range)
    {

        if(from_range.volume() != to_range.volume()){
            throw qk::exception("qk::data::copy_to : Array mismatch.");
        }

        qk::indexer from_idx = from_chunk.indexer(from_range);
        qk::indexer to_idx = to_chunk.indexer(to_range);

        while(true){

            to_chunk[to_idx] = from_chunk[from_idx];

            to_idx.next();
            from_idx.next();

            if(!to_idx.exists()){
                break;
            }
        }
    }

    inline void
    full_copy(const qk::data::datachunk & from_chunk, qk::data::datachunk & to_chunk)
    {

        if(from_chunk.volume() != to_chunk.volume()){
            throw qk::exception("qk::data::full_copy : Array mismatch.");
        }

        std::memcpy(to_chunk.data(),from_chunk.data(),from_chunk.volume()*sizeof(double));

    }

}
}

#endif // _qk_data_functions_H
