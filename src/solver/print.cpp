#include "print.h"

// QK includes
#include "variable/variable.h"

// STL include
#include <iostream>

namespace qk
{
namespace solver
{

void print::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{

    std::ostream & os = std::cout;

    for(const qk::variable::variable_id & var_id : _input_variable_ids){
        const qk::variable::variable & var = variable_manager.output_variable(var_id);

        os << "Printing Variable '"<<var.name()<<"'\n";

        for (qk::indexer chunk_idx = var.indexer(); chunk_idx.exists(); ++chunk_idx) {
            const qk::data::extended_datachunk & var_chunk = var[chunk_idx];

            for(qk::indexer idx = var_chunk.indexer(var_chunk.internal_range()); idx.exists(); ++idx){

                {

                    os << "Index (";
                    for(int i = 0; i < idx.num_dims(); i++){
                        os << idx[i];
                        if(i!=idx.num_dims()-1){
                            os << ", ";
                        }
                    }
                    os << ") : " << var_chunk[idx] << "\n";
                }
            }

        }
    }

}

}
}
