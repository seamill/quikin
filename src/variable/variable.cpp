#include "variable/variable.h"

namespace qk
{
namespace variable
{

variable::variable():
    _name("void")
{

}

variable::~variable()
{

}

variable::variable(const std::string & name,
    const std::vector<std::string> & component_names,
    const qk::basis::basis & basis,
    const qk::range & global_data_range) :
        _name(name),
        _component_names(component_names)
{
    qk::range data_range = global_data_range;
    data_range.extrude(0, component_names.size());
    this->resize(data_range, basis);
}

}
}

