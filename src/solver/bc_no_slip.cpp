#include "bc_no_slip.h"

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
bc_no_slip::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{

    for(int vi = 0; vi < _output_variable_ids.size(); vi++){
        qk::variable::variable & var = variable_manager.output_variable(_output_variable_ids[vi]);

        for(qk::indexer chunk_idx = var.indexer(); chunk_idx.exists(); chunk_idx.next()){

            qk::data::extended_datachunk & chunk = var[chunk_idx];

            for(const int & i : _dims){

                const qk::range & lower_internal_range = chunk.lower_internal_range(i);
                const qk::range & lower_external_range = chunk.lower_external_range(i);

                const qk::range & upper_internal_range = chunk.upper_internal_range(i);
                const qk::range & upper_external_range = chunk.upper_external_range(i);

                // Lower boundary
                for(int j = 0; j<chunk.num_ghost_layers(); ++j){
                    const int from_slice = lower_internal_range.lower(i)+j;
                    const int to_slice = lower_external_range.upper(i)-1-j;

                    qk::range from_rng(lower_internal_range);
                    from_rng.set(i,from_slice,from_slice+1);

                    qk::range to_rng(lower_external_range);
                    to_rng.set(i,to_slice,to_slice+1);

                    qk::data::copy_to(chunk, from_rng, chunk, to_rng);
                }

                // Upper boundary
                for(int j = 0; j<chunk.num_ghost_layers(); ++j){
                    const int from_slice = upper_internal_range.upper(i)-1-j;
                    const int to_slice = upper_external_range.lower(i)+j;

                    qk::range from_rng(upper_internal_range);
                    from_rng.set(i,from_slice,from_slice+1);

                    qk::range to_rng(upper_external_range);
                    to_rng.set(i,to_slice,to_slice+1);

                    qk::data::copy_to(chunk, from_rng, chunk, to_rng);
                }

                // Reverse the normal velocity field

                // Lower boundary
                {
                    qk::range mom_rng(lower_external_range);
                    mom_rng.set(lower_external_range.num_dims()-1, 1,4);
                    chunk.multiply_by_scalar(mom_rng, -1.);
                }

                // Upper boundary
                {
                    qk::range mom_rng(upper_external_range);
                    mom_rng.set(upper_external_range.num_dims()-1, 1,4);
                    chunk.multiply_by_scalar(mom_rng, -1.);
                }
            }
        }
    }

}

}
}
