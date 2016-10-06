#ifndef _qk_data_datachunk_H
#define _qk_data_datachunk_H

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

class datachunk: public qk::indexer_interface<double>
{
public:

    datachunk();
    virtual ~datachunk();

    datachunk & operator=(const datachunk & chunk);

    virtual void write_VTK(std::ofstream & file, const qk::range & range) const;

    virtual void resize(const qk::range & mesh_range, const qk::range & data_range);

private:

    qk::range _mesh_range;
    qk::range _data_range;

    using qk::indexer_interface<double>::resize;

};

}
}

#endif // _qk_data_datachunk_H
