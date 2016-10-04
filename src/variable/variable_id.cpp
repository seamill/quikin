#include "variable/variable_id.h"

namespace qk
{
namespace variable
{

int variable_id::_tag = 0;

variable_id::variable_id()
{

}

variable_id::~variable_id()
{

}

variable_id::variable_id(const std::string & name,
    const std::vector<std::string> & component_names,
    const qk::basis::basis & basis) :
        _name(name),
        _component_names(component_names),
        _basis(basis)
{

}

}
}

