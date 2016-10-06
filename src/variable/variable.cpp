#include "variable/variable.h"

// STL includes
#include <fstream>

// QK includes
#include "lib/indexer.h"
#include "basis/volume_average.h"
#include "basis/basis.h"
#include "lib/range.h"

namespace qk
{
namespace variable
{

variable::variable() :
        _name("void")
{

}

variable::~variable()
{

}

variable::variable(const std::string & name,
    const std::vector<std::string> & component_names,
    const qk::basis::basis & basis,
    const qk::range & global_mesh_range,
    const int num_ghost_layers) :
        _name(name),
        _component_names(component_names),
        qk::data::dataset(global_mesh_range, basis, component_names.size(), num_ghost_layers)
{

}

void variable::write_vtk(std::ofstream & file) const
{
    if (this->basis() == qk::basis::volume_average()) {

        file << "CELL_DATA " << _mesh_range_global.volume() * _data_range.volume() << "\n";

        for (int i = 0; i < _component_names.size(); i++) {
            const std::string & component_name = _component_names[i];
            file << "SCALARS " << _name << "." << component_name << " float 1\n";
            file << "LOOKUP_TABLE default\n";

            bool single_chunk = true;

            for (qk::indexer dc_idx = indexer(); dc_idx.exists(); dc_idx.next()) {
                const qk::data::extended_datachunk & data = (*this)[dc_idx];
                qk::range internal_range = data.internal_range();
                internal_range.set(internal_range.num_dims() - 1, i, i + 1);

                int num = 0;
                for (qk::indexer idx = data.indexer(internal_range); idx.exists(); idx.next()) {
                    file << double(data[idx]) << "\n";
                    num++;
                }
//                std::cout << "Expected " << _data_range_global.volume() << " elements or "<<internal_range.volume() <<" elements, got " << num << std::endl;

                if (!single_chunk) {
                    throw qk::exception("qk::variable::variable::write_vtk : Writer currently limited to a single chunk.");
                }
                single_chunk = false;

            }
        }

    } else {
        throw qk::exception("qk::variable::variable::write_vtk : Basis not yet supported.");
    }
}

}
}

