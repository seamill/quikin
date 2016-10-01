#include <data/dataset.h>
#include <grid/grid.h>
#include <iostream>

// STL includes
#include <cmath>
#include <string>
#include <ctime>
#include <chrono>
#include <sstream>

// QK include
#include "solvers/qksolver.h"
#include "solvers/qksolver_advection.h"

int advection_test()
{

    // Define solver stuff
    const double velocity[] = {0.25,0.25};
    const int num_frames = 10;
    const int num_steps_per_frame = 10;
    const double time_end = 5.0;
    const double time_dt = time_end / double(num_frames * num_steps_per_frame);
    const std::string work_directory = "/home/smiller/localStorage/kinetic/quikin/data";

    // Find a better dt

    // Define domain
    const int num_dims = 2;
    const int dims[] = {256,256};

    qk::range domain_range(num_dims, dims);

    const double startxs[] = {-0.5,-0.5};
    const double widths[] = {1.0,1.0};

    qk::grid::rectilinear grid(domain_range, startxs, widths);

    qk::solver::advection * p_solver = new qk::solver::advection();
    p_solver->setup(grid, velocity);

    qk::solver::solver & solver = *p_solver;

    double time = 0.0;
    double average_delta_time = 0;

    {
        std::stringstream ss;
        ss << work_directory << "/frame";
        std::cout << "Exporting Frame " << 0 << std::endl;
        solver.write_VTK(ss.str(), "_0.vtk");
    }

    for(int i = 0; i < num_frames; i++){

        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();

        for(int j = 0; j < num_steps_per_frame; j++){
            // Run a dt step
            solver.advance(time, time_dt);
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

            solver.write_VTK(ss.str(), ss2.str());
            std::cout << " (plotting took " << 1.e-6 * double(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - later).count()) << " [s])\n";
        }
    }
    std::cout << "Solver ran with average frame time evaluation of " << average_delta_time << std::endl;

    return EXIT_SUCCESS;
}

int main()
{
    return advection_test();
}
