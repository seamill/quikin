#ifndef _qk_data_shared_datachunk_H
#define _qk_data_shared_datachunk_H

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

class shared_datachunk:
        public datachunk
{
public:

	shared_datachunk();
    ~shared_datachunk();

//    void setup(const QKRange & range);

//    void sendData();
//    void recvData();
//    void recvData_nb();
//    void recvData_b();

protected:

//    int _tag;
//    int _sharedWithRank;

private:


};

}
}

#endif // _qk_data_shared_datachunk_H
