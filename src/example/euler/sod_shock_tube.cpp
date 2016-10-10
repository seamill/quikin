#include "sod_shock_tube.h"




// QK include
#include "basis/volume_average.h"
#include "variable/variable_id.h"
#include "variable/variable_manager.h"
#include "grid/rectilinear.h"
#include "lib/exception.h"

#include "solver/ssprk3.h"
#include "solver/euler.h"
#include "solver/shock_tube.h"
#include "solver/bc_no_slip.h"
#include "solver/bc_periodic.h"
#include "solver/swap.h"
#include "solver/fill.h"
#include "solver/print.h"

// STL includes
#include <cmath>
#include <string>
#include <ctime>
#include <chrono>
#include <sstream>
#include <iostream>

namespace qk
{
namespace example
{
namespace euler
{

void
sod_shock_tube(const std::string & work_directory)
{

    const double gamma = 1.4;

    // Define solver stuff
    const int num_dims = 2;
    const int num_frames = 100;
    const int num_steps_per_frame = 5;
    const double time_end = 0.5;
    const double time_dt = time_end / double(num_frames * num_steps_per_frame);
    std::vector<std::string> component_names;
    component_names.push_back("rho");
    component_names.push_back("px");
    component_names.push_back("py");
    component_names.push_back("pz");
    component_names.push_back("e");
    std::vector<int> no_slip_dims(1,0);
    std::vector<int> periodic_dims;
    for(int i=1;i<num_dims;++i){
      periodic_dims.push_back(i);
    }

    // Find a better dt

    // Define domain
    qk::grid::rectilinear grid;
    {
        const int dims[] = { 500,2,2};
        qk::range domain_range(num_dims, dims);
        const double startxs[] = { -0.5,-.5,-.5};
        const double widths[] = { 1.0,1.,1.};

        grid = qk::grid::rectilinear(domain_range, startxs, widths);
    }

    qk::basis::basis basis = qk::basis::volume_average();

    qk::variable::variable_id var("fluid", component_names, basis);
    qk::variable::variable_id var_1("fluid_1", component_names, basis);
    qk::variable::variable_id var_2("fluid_2", component_names, basis);
    qk::variable::variable_id var_3("fluid_3", component_names, basis);
    qk::variable::variable_id rhs("rhs", component_names, basis);

    qk::solver::ssprk3 ssprk3;
    {
        ssprk3.add_input_variable(var);
        ssprk3.add_input_variable(rhs);
        ssprk3.add_output_variable(var_1);
        ssprk3.add_output_variable(var_2);
        ssprk3.add_output_variable(var_3);
    }

    std::vector<qk::solver::euler> eulers;
    {
        qk::solver::euler euler;
        euler.setup(gamma);
        euler.add_input_variable(var);
        euler.add_output_variable(rhs);
        eulers.push_back(euler);
    }
    {
        qk::solver::euler euler;
        euler.setup(gamma);
        euler.add_input_variable(var_1);
        euler.add_output_variable(rhs);
        eulers.push_back(euler);
    }

    {
        qk::solver::euler euler;
        euler.setup(gamma);
        euler.add_input_variable(var_2);
        euler.add_output_variable(rhs);
        eulers.push_back(euler);
    }

    qk::solver::shock_tube ic;
    {
        ic.add_output_variable(var);
    }

    std::vector<qk::solver::bc_no_slip> no_slip_bcs;
    {
        qk::solver::bc_no_slip bc;
        bc.setup(no_slip_dims);
        bc.add_output_variable(var);
        no_slip_bcs.push_back(bc);
    }
    {
        qk::solver::bc_no_slip bc;
        bc.setup(no_slip_dims);
        bc.add_output_variable(var_1);
        no_slip_bcs.push_back(bc);
    }
    {
        qk::solver::bc_no_slip bc;
        bc.setup(no_slip_dims);
        bc.add_output_variable(var_2);
        no_slip_bcs.push_back(bc);
    }

    std::vector<qk::solver::bc_periodic> periodic_bcs;
    {
        qk::solver::bc_periodic bc;
        bc.setup(periodic_dims);
        bc.add_output_variable(var);
        periodic_bcs.push_back(bc);
    }
    {
        qk::solver::bc_periodic bc;
        bc.setup(periodic_dims);
        bc.add_output_variable(var_1);
        periodic_bcs.push_back(bc);
    }
    {
        qk::solver::bc_periodic bc;
        bc.setup(periodic_dims);
        bc.add_output_variable(var_2);
        periodic_bcs.push_back(bc);
    }

    qk::solver::fill reset_rhs;
    {
        reset_rhs.setup(0.);
        reset_rhs.add_output_variable(rhs);
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
                no_slip_bcs[k].solve(variable_manager);
                if(num_dims > 1){
                    periodic_bcs[k].solve(variable_manager);
                }
                reset_rhs.solve(variable_manager);
                eulers[k].solve(variable_manager);
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
