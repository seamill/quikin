#include "wave_1D.h"


// STL includes
#include <cmath>
#include <string>
#include <ctime>
#include <chrono>
#include <sstream>
#include <iostream>

// QK include
#include "basis/volume_average.h"
#include "variable/variable_id.h"
#include "variable/variable_manager.h"
#include "grid/rectilinear.h"
#include "lib/exception.h"


#include "solver/ssprk3.h"
#include "solver/maxwell.h"
#include "solver/electromagnetic_wave_1D.h"
#include "solver/bc_periodic.h"
#include "solver/swap.h"
#include "solver/fill.h"

namespace qk
{
namespace example
{
namespace maxwell
{

void
wave_1D(const std::string & work_directory)
{

    // Define solver stuff
    const double c = 1;
    const double k = 2*3.14159263;
    const double Bx=1,By=1,Bz=1;
    const int num_frames = 100;
    const int num_steps_per_frame = 10;
    const double time_end = 1.0;
    const double time_dt = time_end / double(num_frames * num_steps_per_frame);
    std::vector<std::string> component_names;
    component_names.push_back("Ex");
    component_names.push_back("Ey");
    component_names.push_back("Ez");
    component_names.push_back("Bx");
    component_names.push_back("By");
    component_names.push_back("Bz");

    const int num_dims = 3;
    std::vector<int> periodic_dims;
    for(int i=0;i<num_dims;i++){periodic_dims.push_back(i);}

    // Find a better dt

    // Define domain
    qk::grid::rectilinear grid;

    if(num_dims==1){
        const int dims[] = { 128};
        const double startxs[] = { -0.5};
        const double widths[] = { 1.0};
        grid = qk::grid::rectilinear(qk::range(num_dims, dims), startxs, widths);
    }
    if(num_dims==2){
        const int dims[] = { 2,128};
        const double startxs[] = { -0.5,-.5};
        const double widths[] = { 1.0,1.};
        grid = qk::grid::rectilinear(qk::range(num_dims, dims), startxs, widths);
    }
    if(num_dims==3){
        const int dims[] = {2,2,128};
        const double startxs[] = {-0.5,-.5,-.5};
        const double widths[] = {1.,1.,1.};
        grid = qk::grid::rectilinear(qk::range(num_dims, dims), startxs, widths);
    }


    qk::basis::basis basis = qk::basis::volume_average();

    qk::variable::variable_id var("field", component_names, basis);
    qk::variable::variable_id var_1("field_1", component_names, basis);
    qk::variable::variable_id var_2("field_2", component_names, basis);
    qk::variable::variable_id var_3("field_3", component_names, basis);
    qk::variable::variable_id rhs("rhs", component_names, basis);

    qk::solver::ssprk3 ssprk3;

    ssprk3.add_input_variable(var);
    ssprk3.add_input_variable(rhs);

    ssprk3.add_output_variable(var_1);
    ssprk3.add_output_variable(var_2);
    ssprk3.add_output_variable(var_3);

    std::vector<qk::solver::maxwell> maxwells;
    {
        qk::solver::maxwell maxwell;
        maxwell.setup(c);
        maxwell.add_input_variable(var);
        maxwell.add_output_variable(rhs);
        maxwells.push_back(maxwell);
    }
    {
        qk::solver::maxwell maxwell;
        maxwell.setup(c);
        maxwell.add_input_variable(var_1);
        maxwell.add_output_variable(rhs);
        maxwells.push_back(maxwell);
    }
    {
        qk::solver::maxwell maxwell;
        maxwell.setup(c);
        maxwell.add_input_variable(var_2);
        maxwell.add_output_variable(rhs);
        maxwells.push_back(maxwell);
    }

    qk::solver::electromagnetic_wave_1D ic;
    {
        ic.add_output_variable(var);
        ic.setup(num_dims-1,c,k,Bx,By,Bz);
    }

    std::vector<qk::solver::bc_periodic> bcs;
    {
        qk::solver::bc_periodic bc;
        bc.setup(periodic_dims);
        bc.add_output_variable(var);
        bcs.push_back(bc);
    }
    {
        qk::solver::bc_periodic bc;
        bc.setup(periodic_dims);
        bc.add_output_variable(var_1);
        bcs.push_back(bc);
    }
    {
        qk::solver::bc_periodic bc;
        bc.setup(periodic_dims);
        bc.add_output_variable(var_2);
        bcs.push_back(bc);
    }

    qk::solver::fill fill;
    {
        fill.setup(0.);
        fill.add_output_variable(rhs);
    }

    qk::solver::swap swap;
    {
        swap.add_output_variable(var);
        swap.add_output_variable(var_3);
    }

    qk::variable::variable_manager variable_manager(grid);

    double average_delta_time = 0;

    ic.solve(variable_manager);

    std::vector<qk::variable::variable_id> write_vars;
    write_vars.push_back(var);

    {
        std::cout << "Exporting Frame 0\n";
        variable_manager.write_vtk(work_directory + "/frame", "_0.vtk", write_vars);
    }

    double time = 0.;

    variable_manager.set_time(time);
    variable_manager.set_dt(time_dt);

    for (int i = 0; i < num_frames; i++) {

        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();

        for (int j = 0; j < num_steps_per_frame; j++) {
            // Run a dt step
            for(int k = 0; k < 3; k++){
                bcs[k].solve(variable_manager);
                fill.solve(variable_manager);
                maxwells[k].solve(variable_manager);
                ssprk3.solve(variable_manager,k);
            }
            swap.solve(variable_manager);
            time += time_dt;
            variable_manager.set_time(time);
        }

        std::chrono::steady_clock::time_point later = std::chrono::steady_clock::now();

        const double delta_time = 1.e-6 * double(std::chrono::duration_cast<std::chrono::microseconds>(later - now).count());
        average_delta_time += delta_time / double(num_frames);

        {
            printf("Exporting Frame %3i after time %1.2e [s]",i+1,delta_time);
            std::stringstream ss, ss2;
            ss << work_directory << "/frame";
            ss2 << "_" << i + 1 << ".vtk";
            variable_manager.write_vtk(ss.str(), ss2.str(), write_vars);
            printf(" (plotting took %1.2e [s])\n",1.e-6 * double(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - later).count()));
        }
    }
    std::cout << "Solver ran with average frame time evaluation of " << average_delta_time << std::endl;

}

}
}
}
