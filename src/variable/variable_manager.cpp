#include "variable/variable_manager.h"

// STL includes
#include <fstream>

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

void
variable_manager::add_variable(const qk::variable::variable_id & var)
{
    if(_variables.find(var) != _variables.end()){
        throw qk::exception("qk::variable::variable_manager::add_variable : Variable already exists.");
    }
    _variables[var] = new qk::variable::variable(var.name(), var.component_names(), var.basis(), _grid);
}

void variable_manager::write_vtk(const std::string & prefix,
    const std::string & suffix,
    const std::vector<qk::variable::variable_id> & variables) const
{

    if(variables.size() > 0){

        const qk::basis::basis & basis = variables[0].basis();
        for(const auto & var_id : variables){
            if(basis != var_id.basis()){
                throw qk::exception("qk::variable::variable_manager::write_VTK : Current limitation - all variables must share a basis.");
            }
        }

        // TODO: Add range information
        std::string rangeInfo = "";

        std::string filename = prefix + rangeInfo + suffix;

        std::ofstream file(filename.c_str());
        if(!file){
            throw qk::exception("qk::solver:advection::write_VTK : File '" + filename + "' could not be created");
        }

        file << "# vtk DataFile Version 2.0\n";
        file << "quikin output\n";
        file << "ASCII\n";

        _grid.write_vtk(file, basis);

        file << "CELL_DATA " << _grid.volume() << "\n";

        // Iterate through the domains in n
        for(const qk::variable::variable_id & var_id : variables){
            const auto var_itr = _variables.find(var_id);
            if(var_itr == _variables.end()){
                throw qk::exception("qk::variable::variable_manager::write_VTK : Variable '"+var_id.name()+"' does not exist.");
            }

            var_itr->second->write_vtk(file);

        }

        file.close();
    }

}

}
}

