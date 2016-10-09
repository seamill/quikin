#include "gem_challenge.h"

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
#include "solver/rk1.h"
#include "solver/euler.h"
#include "solver/maxwell.h"
#include "solver/current_source.h"
#include "solver/lorentz_force.h"
#include "solver/gem_challenge.h"
#include "solver/bc_periodic.h"
#include "solver/bc_no_slip.h"
#include "solver/bc_conducting_wall.h"
#include "solver/swap.h"
#include "solver/fill.h"
#include "solver/print.h"

namespace qk
{
namespace example
{
namespace plasma
{

void
gem_challenge()
{

    typedef qk::solver::ssprk3 time_integrator;


    const double e = 1.60217662e-19;
    const double me = 9.10938356e-31;
    const double mp = 1.6726219e-27;
    const double eps0 = 8.85418782e-12;
    const double mu0 = 1.25663706e-6;
    const double c = 299792458.0;
    const double gamma = 5./3.;
    const double pi = 3.14159263;



    const double di = 1.;
    const double qe = -e;
    const double qi = e;

    const double mi_me = 25.;
    const double lam_di = 0.5;
    const double n1_n0 = 0.2;
    const double Te_Ti = 0.2;
    const double vse_c = 0.1;
    const double Lx = 4 * pi * di;

//    const double Lx = 2 * pi * di;
//    const double mi_me = 1.;
//    const double Te_Ti = 1.0;

    const double mi = mi_me * me;
    const double n0 = c*c*mi*eps0/(e*e*di*di);
    const double n1 = n1_n0 * n0;
    const double lam = lam_di * di;

    const double vse = vse_c * c;
    const double Te = vse * vse * me / gamma;
    const double Ti = Te / Te_Ti;
    const double P0 = n0 * (Te + Ti);
    const double B0 = std::sqrt(2*mu0*P0);

    const double ue0 = 2*Te/qe/lam/B0;
    const double ui0 = 2*Ti/qi/lam/B0;

    const double wce = e*B0/me;
    const double wci = e*B0/mi;

    const double wpe = std::sqrt(n0*e*e/eps0/me);
    const double wpi = std::sqrt(n0*e*e/eps0/mi);

    const int Nx = 128;

    std::vector<double> masses;
    masses.push_back(me);
    masses.push_back(mi);
    std::vector<double> charges;
    charges.push_back(-e);
    charges.push_back( e);

    // Define solver stuff
    const int num_frames = 1000;
    const int num_steps_per_frame = 20;
    const double time_end = 40./wci;
    const double time_dt = time_end / double(num_frames * num_steps_per_frame);
    const std::vector<int> periodic_dims(1,1);
    const std::vector<int> walls_dims(1,0);
    const std::string work_directory = "/Users/seamill/local_storage/quikin/data";
    std::vector<std::string> fluid_component_names;
    fluid_component_names.push_back("rho");
    fluid_component_names.push_back("px");
    fluid_component_names.push_back("py");
    fluid_component_names.push_back("pz");
    fluid_component_names.push_back("e");

    std::vector<std::string> field_component_names;
    field_component_names.push_back("Ex");
    field_component_names.push_back("Ey");
    field_component_names.push_back("Ez");
    field_component_names.push_back("Bx");
    field_component_names.push_back("By");
    field_component_names.push_back("Bz");

    // Define domain
    qk::grid::rectilinear grid;
    {
        const int num_dims = 2;
        const int dims[] = {Nx,2*Nx};
        qk::range domain_range(num_dims, dims);
        const double startxs[] = {-Lx/2.,-Lx};
        const double widths[] = {Lx,2*Lx};

        grid = qk::grid::rectilinear(domain_range, startxs, widths);
    }




    const double dx = grid.dx(0);
    std::cout << "CFL Speed of light: " << time_dt * c / dx << std::endl;
    std::cout << "CFL Speed of sound (e): " << time_dt * vse / dx << std::endl;
    std::cout << "dt * wpe: " << time_dt * wpe << std::endl;
    std::cout << "dt * wce: " << time_dt * wce << std::endl;
    std::cout << "dt * wpi: " << time_dt * wpi << std::endl;
    std::cout << "dt * wci: " << time_dt * wci << std::endl;




    qk::basis::basis basis = qk::basis::volume_average();

    std::vector<qk::variable::variable_id> electrons;
    {
        electrons.push_back(qk::variable::variable_id("electrons"    , fluid_component_names, basis));
        electrons.push_back(qk::variable::variable_id("electrons_1"  , fluid_component_names, basis));
        electrons.push_back(qk::variable::variable_id("electrons_2"  , fluid_component_names, basis));
        electrons.push_back(qk::variable::variable_id("electrons_3"  , fluid_component_names, basis));
        electrons.push_back(qk::variable::variable_id("electrons_rhs", fluid_component_names, basis));
    }

    std::vector<qk::variable::variable_id> ions;
    {
        ions.push_back(qk::variable::variable_id("ions"    , fluid_component_names, basis));
        ions.push_back(qk::variable::variable_id("ions_1"  , fluid_component_names, basis));
        ions.push_back(qk::variable::variable_id("ions_2"  , fluid_component_names, basis));
        ions.push_back(qk::variable::variable_id("ions_3"  , fluid_component_names, basis));
        ions.push_back(qk::variable::variable_id("ions_rhs", fluid_component_names, basis));
    }

    std::vector<qk::variable::variable_id> fields;
    {
        fields.push_back(qk::variable::variable_id("field"    , field_component_names, basis));
        fields.push_back(qk::variable::variable_id("field_1"  , field_component_names, basis));
        fields.push_back(qk::variable::variable_id("field_2"  , field_component_names, basis));
        fields.push_back(qk::variable::variable_id("field_3"  , field_component_names, basis));
        fields.push_back(qk::variable::variable_id("field_rhs", field_component_names, basis));
    }

    time_integrator electrons_ti;
    {
        electrons_ti.add_input_variable(electrons[0]);
        electrons_ti.add_input_variable(electrons[4]);
        electrons_ti.add_output_variable(electrons[1]);
        electrons_ti.add_output_variable(electrons[2]);
        electrons_ti.add_output_variable(electrons[3]);
    }

    time_integrator ions_ti;
    {
        ions_ti.add_input_variable(ions[0]);
        ions_ti.add_input_variable(ions[4]);
        ions_ti.add_output_variable(ions[1]);
        ions_ti.add_output_variable(ions[2]);
        ions_ti.add_output_variable(ions[3]);
    }

    time_integrator field_ti;
    {
        field_ti.add_input_variable(fields[0]);
        field_ti.add_input_variable(fields[4]);
        field_ti.add_output_variable(fields[1]);
        field_ti.add_output_variable(fields[2]);
        field_ti.add_output_variable(fields[3]);
    }


    std::vector<qk::solver::euler> electron_eulers;
    {
        qk::solver::euler ss;
        ss.setup(gamma);
        ss.add_input_variable(electrons[0]);
        ss.add_output_variable(electrons[4]);
        electron_eulers.push_back(ss);
    }
    {
        qk::solver::euler ss;
        ss.setup(gamma);
        ss.add_input_variable(electrons[1]);
        ss.add_output_variable(electrons[4]);
        electron_eulers.push_back(ss);
    }
    {
        qk::solver::euler ss;
        ss.setup(gamma);
        ss.add_input_variable(electrons[2]);
        ss.add_output_variable(electrons[4]);
        electron_eulers.push_back(ss);
    }

    std::vector<qk::solver::lorentz_force> electron_lorentz_forces;
    {
        qk::solver::lorentz_force ss;
        ss.setup(qe,me);
        ss.add_input_variable(electrons[0]);
        ss.add_input_variable(fields[0]);
        ss.add_output_variable(electrons[4]);
        electron_lorentz_forces.push_back(ss);
    }
    {
        qk::solver::lorentz_force ss;
        ss.setup(qe,me);
        ss.add_input_variable(electrons[1]);
        ss.add_input_variable(fields[1]);
        ss.add_output_variable(electrons[4]);
        electron_lorentz_forces.push_back(ss);
    }
    {
        qk::solver::lorentz_force ss;
        ss.setup(qe,me);
        ss.add_input_variable(electrons[2]);
        ss.add_input_variable(fields[2]);
        ss.add_output_variable(electrons[4]);
        electron_lorentz_forces.push_back(ss);
    }

    std::vector<qk::solver::euler> ion_eulers;
    {
        qk::solver::euler ss;
        ss.setup(gamma);
        ss.add_input_variable(ions[0]);
        ss.add_output_variable(ions[4]);
        ion_eulers.push_back(ss);
    }
    {
        qk::solver::euler ss;
        ss.setup(gamma);
        ss.add_input_variable(ions[1]);
        ss.add_output_variable(ions[4]);
        ion_eulers.push_back(ss);
    }
    {
        qk::solver::euler ss;
        ss.setup(gamma);
        ss.add_input_variable(ions[2]);
        ss.add_output_variable(ions[4]);
        ion_eulers.push_back(ss);
    }

    std::vector<qk::solver::lorentz_force> ion_lorentz_forces;
    {
        qk::solver::lorentz_force ss;
        ss.setup(qi,mi);
        ss.add_input_variable(ions[0]);
        ss.add_input_variable(fields[0]);
        ss.add_output_variable(ions[4]);
        ion_lorentz_forces.push_back(ss);
    }
    {
        qk::solver::lorentz_force ss;
        ss.setup(qi,mi);
        ss.add_input_variable(ions[1]);
        ss.add_input_variable(fields[1]);
        ss.add_output_variable(ions[4]);
        ion_lorentz_forces.push_back(ss);
    }
    {
        qk::solver::lorentz_force ss;
        ss.setup(qi,mi);
        ss.add_input_variable(ions[2]);
        ss.add_input_variable(fields[2]);
        ss.add_output_variable(ions[4]);
        ion_lorentz_forces.push_back(ss);
    }


    std::vector<qk::solver::maxwell> field_maxwells;
    {
        qk::solver::maxwell ss;
        ss.setup(c);
        ss.add_input_variable(fields[0]);
        ss.add_output_variable(fields[4]);
        field_maxwells.push_back(ss);
    }
    {
        qk::solver::maxwell ss;
        ss.setup(c);
        ss.add_input_variable(fields[1]);
        ss.add_output_variable(fields[4]);
        field_maxwells.push_back(ss);
    }
    {
        qk::solver::maxwell ss;
        ss.setup(c);
        ss.add_input_variable(fields[2]);
        ss.add_output_variable(fields[4]);
        field_maxwells.push_back(ss);
    }

    std::vector<qk::solver::current_source> field_current_sources;
    {
        qk::solver::current_source ss;
        ss.setup(eps0, charges, masses);
        ss.add_input_variable(electrons[0]);
        ss.add_input_variable(ions[0]);
        ss.add_output_variable(fields[4]);
        field_current_sources.push_back(ss);
    }
    {
        qk::solver::current_source ss;
        ss.setup(eps0, charges, masses);
        ss.add_input_variable(electrons[1]);
        ss.add_input_variable(ions[1]);
        ss.add_output_variable(fields[4]);
        field_current_sources.push_back(ss);
    }
    {
        qk::solver::current_source ss;
        ss.setup(eps0, charges, masses);
        ss.add_input_variable(electrons[2]);
        ss.add_input_variable(ions[2]);
        ss.add_output_variable(fields[4]);
        field_current_sources.push_back(ss);
    }

    qk::solver::gem_challenge ic;
    {
        ic.add_output_variable(electrons[0]);
        ic.add_output_variable(ions[0]);
        ic.add_output_variable(fields[0]);
        ic.setup(lam,Lx,n0,n1,gamma,qe,me,ue0,Te,qi,mi,ui0,Ti,B0);
    }

    std::vector<qk::solver::bc_no_slip> fluid_wall_bcs;
    {
        qk::solver::bc_no_slip bc;
        bc.setup(walls_dims);
        bc.add_output_variable(electrons[0]);
        bc.add_output_variable(ions[0]);
        fluid_wall_bcs.push_back(bc);
    }
    {
        qk::solver::bc_no_slip bc;
        bc.setup(walls_dims);
        bc.add_output_variable(electrons[1]);
        bc.add_output_variable(ions[1]);
        fluid_wall_bcs.push_back(bc);
    }
    {
        qk::solver::bc_no_slip bc;
        bc.setup(walls_dims);
        bc.add_output_variable(electrons[2]);
        bc.add_output_variable(ions[2]);
        fluid_wall_bcs.push_back(bc);
    }

    std::vector<qk::solver::bc_conducting_wall> field_bcs;
    {
        qk::solver::bc_conducting_wall bc;
        bc.setup(walls_dims);
        bc.add_output_variable(fields[0]);
        field_bcs.push_back(bc);
    }
    {
        qk::solver::bc_conducting_wall bc;
        bc.setup(walls_dims);
        bc.add_output_variable(fields[1]);
        field_bcs.push_back(bc);
    }
    {
        qk::solver::bc_conducting_wall bc;
        bc.setup(walls_dims);
        bc.add_output_variable(fields[2]);
        field_bcs.push_back(bc);
    }

    std::vector<qk::solver::bc_periodic> periodic_bcs;
    {
        qk::solver::bc_periodic bc;
        bc.setup(periodic_dims);
        bc.add_output_variable(electrons[0]);
        bc.add_output_variable(ions[0]);
        bc.add_output_variable(fields[0]);
        periodic_bcs.push_back(bc);
    }
    {
        qk::solver::bc_periodic bc;
        bc.setup(periodic_dims);
        bc.add_output_variable(electrons[1]);
        bc.add_output_variable(ions[1]);
        bc.add_output_variable(fields[1]);
        periodic_bcs.push_back(bc);
    }
    {
        qk::solver::bc_periodic bc;
        bc.setup(periodic_dims);
        bc.add_output_variable(electrons[2]);
        bc.add_output_variable(ions[2]);
        bc.add_output_variable(fields[2]);
        periodic_bcs.push_back(bc);
    }

    qk::solver::fill reset_rhs;
    {
        reset_rhs.setup(0.);
        reset_rhs.add_output_variable(electrons[4]);
        reset_rhs.add_output_variable(ions[4]);
        reset_rhs.add_output_variable(fields[4]);
    }

    qk::solver::swap swap;
    {
        swap.add_output_variable(electrons[0]);
        swap.add_output_variable(electrons[electrons_ti.num_stages()]);
        swap.add_output_variable(ions[0]);
        swap.add_output_variable(ions[ions_ti.num_stages()]);
        swap.add_output_variable(fields[0]);
        swap.add_output_variable(fields[field_ti.num_stages()]);
    }

    std::vector<qk::variable::variable_id> write_vars;
    write_vars.push_back(electrons[0]);
    write_vars.push_back(ions[0]);
    write_vars.push_back(fields[0]);

    qk::variable::variable_manager variable_manager(grid);

    double average_delta_time = 0;

    ic.solve(variable_manager);

    {
        std::cout << "Exporting Frame   0\n";
        variable_manager.write_vtk(work_directory + "/frame", "_0.vtk", write_vars);
    }

    double time = 0.;

    variable_manager.set_time(time);
    variable_manager.set_dt(time_dt);

    for (int i = 0; i < num_frames; i++) {

        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();

        for (int j = 0; j < num_steps_per_frame; j++) {
            // Run a dt step
            for(int k = 0; k < electrons_ti.num_stages(); k++){
                fluid_wall_bcs[k].solve(variable_manager);
                field_bcs[k].solve(variable_manager);
                periodic_bcs[k].solve(variable_manager);
                reset_rhs.solve(variable_manager);
                electron_eulers[k].solve(variable_manager);
                electron_lorentz_forces[k].solve(variable_manager);
                ion_eulers[k].solve(variable_manager);
                ion_lorentz_forces[k].solve(variable_manager);
                field_maxwells[k].solve(variable_manager);
                field_current_sources[k].solve(variable_manager);

                electrons_ti.solve(variable_manager,k);
                ions_ti.solve(variable_manager,k);
                field_ti.solve(variable_manager,k);
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
