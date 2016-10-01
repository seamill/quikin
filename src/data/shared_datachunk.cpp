#include "shared_datachunk.h"

// STL includes
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>

namespace qk
{
namespace data
{

shared_datachunk::shared_datachunk():
    datachunk()
{

}

shared_datachunk::~shared_datachunk()
{
    // Nothing to do - data array is cleared by indexer interface
}


}
}

