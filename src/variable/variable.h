#ifndef _qk_variable_variable_H
#define _qk_variable_variable_H

// STL includes
#include <string>
#include <vector>

// QK includes
#include "basis/basis.h"
#include "data/dataset.h"
#include "lib/range.h"

namespace qk
{
namespace variable
{

class variable: public qk::data::dataset
{
public:

    variable();

    variable(const std::string & name, const std::vector<std::string> & component_names, const qk::basis::basis & basis, const qk::range & global_data_range);

    ~variable();

protected:

    std::string _name;
    std::vector<std::string> _component_names;

};

}
}

#endif // _qk_variable_variable_H
