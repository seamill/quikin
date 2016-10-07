
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
#include "grid/grid.h"
#include "lib/exception.h"
#include "solver/ssprk2_advection.h"
#include "solver/gaussian.h"
#include "solver/bc_periodic.h"
#include "solver/swap.h"

namespace qk
{
namespace example
{
namespace advection
{

void
advect_gaussian()
{

    // Define solver stuff
    std::vector<double> velocity(2, 1.0);
    velocity[0] = -1.;
    velocity[1] = -1.;
    const int num_frames = 100;
    const int num_steps_per_frame = 10;
    const double time_end = 1.0;
    const double time_dt = time_end / double(num_frames * num_steps_per_frame);
    const std::string work_directory = "/Users/seamill/local_storage/quikin/data";

    // Find a better dt

    // Define domain
    qk::grid::rectilinear grid;
    {
        const int num_dims = 2;
        const int dims[] = { 128, 128 };
        qk::range domain_range(num_dims, dims);
        const double startxs[] = { -0.5, -0.5 };
        const double widths[] = { 1.0, 1.0 };

        grid = qk::grid::rectilinear(domain_range, startxs, widths);
    }

    qk::basis::basis basis = qk::basis::volume_average();
//
//    qk::variable::variable density;
//    {
//        std::vector<std::string> component_names("density",1);
//        density = qk::data::variable(component_names, basis, domain_range);
//    }
//
//    qk::solver::advection * p_solver = new qk::solver::advection();
//    p_solver->setup(grid, velocity);

//    qk::solver::solver & solver = *p_solver;

    qk::variable::variable_id var;
    {
        std::vector<std::string> component_names;
        component_names.push_back("rho");
        var = qk::variable::variable_id("variable", component_names, basis);
    }

    qk::variable::variable_id partial_var;
    {
        std::vector<std::string> component_names;
        component_names.push_back("rho");
        partial_var = qk::variable::variable_id("partial_variable", component_names, basis);
    }

    qk::variable::variable_id new_var;
    {
        std::vector<std::string> component_names;
        component_names.push_back("rho");
        new_var = qk::variable::variable_id("new_variable", component_names, basis);
    }

    qk::solver::ssprk2_advection solver;
    {

        std::vector<qk::variable::variable_id> in_vars;
        in_vars.push_back(var);

        std::vector<qk::variable::variable_id> out_vars;
        out_vars.push_back(partial_var);
        out_vars.push_back(new_var);

        solver.setup(velocity);
        solver.set_input_variables(in_vars);
        solver.set_output_variables(out_vars);
    }

    qk::solver::gaussian gaussian;
    {
        std::vector<qk::variable::variable_id> in_vars;
        in_vars.push_back(var);

        std::vector<double> average(3, 0);
        gaussian.set_output_variables(in_vars);
        gaussian.setup(1.0, average, 0.1);
    }

    qk::solver::bc_periodic bc_0;
    {
        std::vector<qk::variable::variable_id> out_vars;
        out_vars.push_back(var);

        bc_0.set_output_variables(out_vars);
    }

    qk::solver::bc_periodic bc_1;
    {
        std::vector<qk::variable::variable_id> out_vars;
        out_vars.push_back(partial_var);

        bc_1.set_output_variables(out_vars);
    }

    qk::solver::swap swap;
    {
        std::vector<qk::variable::variable_id> swap_vars;
        swap_vars.push_back(var);
        swap_vars.push_back(new_var);

        swap.set_output_variables(swap_vars);
    }

    qk::variable::variable_manager variable_manager(grid);

    double average_delta_time = 0;

    gaussian.solve(variable_manager);

    std::vector<qk::variable::variable_id> write_vars;
    write_vars.push_back(var);

    {
        std::cout << "Exporting Frame 0\n";
        variable_manager.write_vtk(work_directory + "/frame", "_0.vtk", write_vars);
        std::cout << "Exporting complete\n";
    }

    double time = 0.;

    variable_manager.set_time(time);
    variable_manager.set_dt(time_dt);

    for (int i = 0; i < num_frames; i++) {

        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();

        for (int j = 0; j < num_steps_per_frame; j++) {
            // Run a dt step
            bc_0.solve(variable_manager);
            solver.solve(variable_manager,qk::solver::ssprk2_advection::STAGE_0);
            bc_1.solve(variable_manager);
            solver.solve(variable_manager,qk::solver::ssprk2_advection::STAGE_1);
            swap.solve(variable_manager);
            time += time_dt;
            variable_manager.set_time(time);
        }

        std::chrono::steady_clock::time_point later = std::chrono::steady_clock::now();

        const double delta_time = 1.e-6 * double(std::chrono::duration_cast<std::chrono::microseconds>(later - now).count());
        average_delta_time += delta_time / double(num_frames);

        {
            std::cout << "Exporting Frame " << i + 1 << " after time " << delta_time << " [s]";
            std::stringstream ss, ss2;
            ss << work_directory << "/frame";
            ss2 << "_" << i + 1 << ".vtk";
            variable_manager.write_vtk(ss.str(), ss2.str(), write_vars);
            std::cout << " (plotting took " << 1.e-6 * double(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - later).count()) << " [s])\n";
        }
    }
    std::cout << "Solver ran with average frame time evaluation of " << average_delta_time << std::endl;

}

}
}
}
