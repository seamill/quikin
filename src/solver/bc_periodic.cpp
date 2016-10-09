#include "bc_periodic.h"

// QK includes
#include "grid/grid.h"
#include "data/functions.h"
#include "variable/variable.h"

// STL include
#include <cmath>
#include <iomanip>

namespace qk
{
namespace solver
{

void
bc_periodic::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{

    for(int vi = 0; vi < _output_variable_ids.size(); vi++){
        qk::variable::variable & var = variable_manager.output_variable(_output_variable_ids[vi]);

        for(qk::indexer chunk_idx = var.indexer(); chunk_idx.exists(); chunk_idx.next()){

            //TODO: This will have to be modified in the future for multiple chunks
            qk::data::extended_datachunk & chunk = var[chunk_idx];

            for(const int & i : _dims){

                qk::data::copy_to(chunk, chunk.lower_internal_range(i), chunk, chunk.upper_external_range(i));
                qk::data::copy_to(chunk, chunk.upper_internal_range(i), chunk, chunk.lower_external_range(i));

            }
        }
    }

}

}
}
