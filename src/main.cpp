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
#include "solver/swap.h"

int advection_test()
{

    // Define solver stuff
    const std::vector<double> velocity(0.25,2);
    const int num_frames = 10;
    const int num_steps_per_frame = 10;
    const double time_end = 5.0;
    const double time_dt = time_end / double(num_frames * num_steps_per_frame);
    const std::string work_directory = "/home/smiller/local_storage/quikin/data";

    // Find a better dt

    // Define domain
    qk::grid::rectilinear grid;
    {
        const int num_dims = 2;
        const int dims[] = {256,256};
        qk::range domain_range(num_dims, dims);
        const double startxs[] = {-0.5,-0.5};
        const double widths[] = {1.0,1.0};

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
        var = qk::variable::variable_id("variable",component_names, basis);
    }

    std::vector<qk::variable::variable_id> vars;
    vars.push_back(var);

    qk::variable::variable_id new_var;
    {
        std::vector<std::string> component_names;
        component_names.push_back("rho");
        var = qk::variable::variable_id("new_variable",component_names, basis);
    }

    qk::solver::ssprk2_advection solver;
    {

        std::vector<qk::variable::variable_id> new_vars;
        new_vars.push_back(new_var);

        solver.setup(velocity);
        solver.set_input_variables(vars);
        solver.set_output_variables(new_vars);
    }

    qk::solver::gaussian gaussian;
    {
        std::vector<double> average(0.,3);
        gaussian.set_output_variables(vars);
        gaussian.setup(1.0,average, 0.1);
    }

    qk::solver::swap swap;
    {
        std::vector<qk::variable::variable_id> swap_vars;
        swap_vars.push_back(var);
        swap_vars.push_back(new_var);

        swap.set_output_variables(swap_vars);
    }

    qk::variable::variable_manager variable_manager(grid);

    double time = 0.0;
    double average_delta_time = 0;

    gaussian.solve(time, variable_manager);

    {
        std::cout << "Exporting Frame 0\n";
        variable_manager.write_VTK(work_directory + "/frame", "_0.vtk", vars);
    }

    for(int i = 0; i < num_frames; i++){

        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();

        for(int j = 0; j < num_steps_per_frame; j++){
            // Run a dt step
            solver.solve(time, variable_manager);
            swap.solve(time, variable_manager);
            time+=time_dt;
        }

        std::chrono::steady_clock::time_point later = std::chrono::steady_clock::now();

        const double delta_time = 1.e-6 * double(std::chrono::duration_cast<std::chrono::microseconds>(later - now).count());
        average_delta_time += delta_time / double(num_frames);

        {
            std::cout << "Exporting Frame " << i+1 << " after time " << delta_time << " [s]";
            std::stringstream ss, ss2;
            ss << work_directory << "/frame";
            ss2 << "_" << i+1 << ".vtk";
            variable_manager.write_VTK(ss.str(), ss2.str(), vars);
            std::cout << " (plotting took " << 1.e-6 * double(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - later).count()) << " [s])\n";
        }
    }
    std::cout << "Solver ran with average frame time evaluation of " << average_delta_time << std::endl;

    return EXIT_SUCCESS;
}

int main()
{
    try{
        advection_test();
    } catch (qk::exception & qke){
        std::cout << "*** quikin exception caught ***\n" << qke.what();
        exit(EXIT_FAILURE);
    }
}
