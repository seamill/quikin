#include "swap.h"

// STL include
#include <cmath>

// QK includes
#include "grid/grid.h"
#include "variable/variable.h"

namespace qk
{
namespace solver
{

swap::swap()
{

}

swap::~swap()
{

}

void
swap::solve(const double time, qk::variable::variable_manager & variable_manager) const
{

    if(_output_variable_ids.size() % 2 != 0){
        throw qk::exception("qk::solver::swap::solve : Must have even number of output variables.");
    }

    for(int i = 0; i < _output_variable_ids.size(); i++){

        qk::variable::variable & a = variable_manager.output_variable(_output_variable_ids[2*i]);
        qk::variable::variable & b = variable_manager.output_variable(_output_variable_ids[2*i+1]);

        a.swap(b);

    }

}

}
}
