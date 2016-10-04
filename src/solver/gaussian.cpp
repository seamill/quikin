#include "gaussian.h"

// STL include
#include <cmath>

// QK includes
#include "grid/grid.h"
#include "variable/variable.h"

namespace qk
{
namespace solver
{

gaussian::gaussian()
{

}

gaussian::~gaussian()
{

}

void
gaussian::setup(const double amplitude, const std::vector<double> & average, const double standard_deviation)
{
    _amplitude = amplitude;
    _average = average;
    _average.resize(3,0.);
    _standard_deviation = standard_deviation;

}

void
gaussian::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{

    const qk::grid::grid & grid = variable_manager.grid();

    double x[3];
    double w[3];

    for(int i = 0; i < _output_variable_ids.size(); i++){
        qk::variable::variable & var = variable_manager.output_variable(_output_variable_ids[i]);

        for(qk::indexer data_idx = var.indexer(); data_idx.exists(); data_idx.next()){
            qk::data::datachunk & data = var[data_idx];
            for(qk::indexer idx = data.indexer(); idx.exists(); idx.next()){
                grid.xc(idx,x);

                w[0] = (x[0] - _average[0]) / _standard_deviation;
                w[1] = (x[1] - _average[1]) / _standard_deviation;
                w[2] = (x[2] - _average[2]) / _standard_deviation;

                data[idx] = _amplitude * std::exp(-w[0]*w[0]-w[1]*w[1]-w[2]*w[2]);
            }

        }

    }

}

}
}
