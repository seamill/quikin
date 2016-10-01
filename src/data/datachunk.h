#ifndef _qk_data_datachunk_H
#define _qk_data_datachunk_H

// STL includes
#include <string>

// QK includes
#include "lib/qkrange.h"
#include "lib/qkindexer.h"
#include "lib/qkindexerinterface.h"

namespace qk
{
namespace data
{

class datachunk:
        public qk::indexer_interface<double>
{
public:

    datachunk();
    virtual ~datachunk();

    datachunk & operator=(const datachunk & dataset);

    void write_VTK(std::ofstream & file, const qk::range & range) const;

};

}
}

#endif // _qk_data_datachunk_H