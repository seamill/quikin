#include "bc_periodic.h"

// STL include
#include <cmath>

// QK includes
#include "grid/grid.h"
#include "data/functions.h"
#include "variable/variable.h"

namespace qk
{
namespace solver
{

bc_periodic::bc_periodic()
{

}

bc_periodic::~bc_periodic()
{

}


void
bc_periodic::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{

    const int num_dims = variable_manager.grid().range().num_dims();

    for(int i = 0; i < _output_variable_ids.size(); i++){
        qk::variable::variable & var = variable_manager.output_variable(_output_variable_ids[i]);

        for(qk::indexer chunk_idx = var.indexer(); chunk_idx.exists(); chunk_idx.next()){

            //TODO: This will have to be modified in the future for multiple chunks
            const qk::data::extended_datachunk & from_chunk = var[chunk_idx];
            qk::data::extended_datachunk & to_chunk = var[chunk_idx];

            std::cout << "Full range = " << from_chunk.range() << std::endl;

            for(int i=0;i<num_dims;i++){

                qk::data::copy_to(from_chunk, from_chunk.lower_internal_range(i), to_chunk, to_chunk.upper_external_range(i));
                qk::data::copy_to(from_chunk, from_chunk.upper_internal_range(i), to_chunk, to_chunk.lower_external_range(i));

                std::cout << "low to up " << i << " = " <<  from_chunk.lower_internal_range(i) << ", " << to_chunk.upper_external_range(i) << std::endl;
                std::cout << "up to low " << i << " = " <<  from_chunk.upper_internal_range(i) << ", " << to_chunk.lower_external_range(i) << std::endl;


            }
        }
    }

}

}
}
