#include "fill.h"

// STL include
#include <cmath>

// QK includes
#include "variable/variable.h"

namespace qk
{
namespace solver
{

void fill::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{

    for(int i = 0; i < _output_variable_ids.size(); i++){
        qk::variable::variable & q = variable_manager.output_variable(_output_variable_ids[i]);

        for (qk::indexer idx = q.indexer(); idx.exists(); ++idx) {
            q[idx].fill(_value);
        }
    }

}

}
}
