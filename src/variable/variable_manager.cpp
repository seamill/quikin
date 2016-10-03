#include "variable/variable_manager.h"

// STL includes

// QK includes
#include "variable/variable.h"

namespace qk
{
namespace variable
{

variable_manager::variable_manager(const qk::grid::grid & grid) :
        _grid(grid)
{

}

variable_manager::~variable_manager()
{
    for (auto & var_pair : _variables) {
        delete var_pair.second;
        var_pair.second = NULL;
    }
    _variables.clear();
}

qk::variable::variable &
variable_manager::output_variable(const qk::variable::variable_id & var)
{
    auto itr = _variables.find(var);
    if (itr == _variables.end()) {
       add_variable(var);
       itr = _variables.find(var);
       if (itr == _variables.end()) {
           qk::exception("qk::variable::variable_manager::output_variable : Couldn't allocate variable.");
       }
    }
    return *(itr->second);
}

const qk::variable::variable &
variable_manager::input_variable(const qk::variable::variable_id & var) const
{
    const auto & itr = _variables.find(var);
    if (itr == _variables.end()) {
       qk::exception("qk::variable::variable_manager::input_variable : Variable '"+var.name()+"' does not exist.");
    }
    return *(itr->second);
}

void variable_manager::add_variable(const qk::variable::variable_id & var)
{
    qk::variable::variable * variable = new qk::variable::variable(var.name(), var.component_names(), var.basis(), _grid);
}

void variable_manager::write_VTK(const std::string & prefix,
    const std::string & suffix,
    const std::vector<qk::variable::variable_id> & variables) const
{

    throw qk::exception("qk::variable::variable_manager::write_VTK : Not yet setup.");

//    // Iterate through the domains in n
//    for(qk::indexer indexer = _n.indexer(); indexer.exists(); indexer.next()){
//        const qk::data::extended_datachunk & chunk = _n[indexer];
//        const qk::range & range = chunk.internal_range();
//
//        // Generate Filename
//
//        // TODO: Add range information
//        std::string rangeInfo = "";
//
//        std::string filename = prefix + rangeInfo + suffix;
//
//        std::ofstream file(filename.c_str());
//        if(!file){
//            throw qk::exception("qk::solver:advection::write_VTK : File " + filename + " could not be created");
//        }
//
//        file << "# vtk DataFile Version 2.0\n";
//        file << "Stuff for dataset\n";
//        file << "ASCII\n";
//
//        //_grid.writeVTK(file, range);
//        chunk.write_VTK(file, range);
//
//    }

}

}
}

