#include "two_fluid_brio_wu.h"

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

#include "solver/ssprk3.h"
#include "solver/euler.h"
#include "solver/maxwell.h"
#include "solver/current_source.h"
#include "solver/lorentz_force.h"
#include "solver/brio_wu_shock_tube.h"
#include "solver/bc_periodic.h"
#include "solver/swap.h"
#include "solver/fill.h"

namespace qk
{
namespace example
{
namespace plasma
{

void
two_fluid_brio_wu()
{

    const double c = 10.;
    const double eps0 = 1.;
    const double gamma = 2.;
    const double qe = -1.0;
    const double qi = 1.0;
    const double me = 0.04;
    const double mi = 1.0;
    const double Pl = 1.0;
    const double Pr = 0.1;
    const double nl = 1.0;
    const double nr = 0.125;
    const double Bx = 0.75;
    const double By = 1.;



    std::vector<double> masses;
    masses.push_back(me);
    masses.push_back(mi);
    std::vector<double> charges;
    charges.push_back(qe);
    charges.push_back(qi);




    // Define solver stuff
    const int num_frames = 100;
    const int num_steps_per_frame = 10;
    const double time_end = 1.;
    const double time_dt = time_end / double(num_frames * num_steps_per_frame);
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


    // Find a better dt

    // Define domain
    qk::grid::rectilinear grid;
    {
        const int num_dims = 2;
        const int dims[] = { 100, 2 };
        qk::range domain_range(num_dims, dims);
        const double startxs[] = { -0.5, -0.5 };
        const double widths[] = { 1.0, 1.0 };

        grid = qk::grid::rectilinear(domain_range, startxs, widths);
    }

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

    qk::solver::ssprk3 electrons_ti;
    {
        electrons_ti.add_input_variable(electrons[0]);
        electrons_ti.add_input_variable(electrons[4]);
        electrons_ti.add_output_variable(electrons[1]);
        electrons_ti.add_output_variable(electrons[2]);
        electrons_ti.add_output_variable(electrons[3]);
    }

    qk::solver::ssprk3 ions_ti;
    {
        ions_ti.add_input_variable(ions[0]);
        ions_ti.add_input_variable(ions[4]);
        ions_ti.add_output_variable(ions[1]);
        ions_ti.add_output_variable(ions[2]);
        ions_ti.add_output_variable(ions[3]);
    }

    qk::solver::ssprk3 field_ti;
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


    qk::solver::brio_wu_shock_tube ic;
    {
        ic.add_output_variable(electrons[0]);
        ic.add_output_variable(ions[0]);
        ic.add_output_variable(fields[0]);
        ic.setup(gamma,me,mi,nl,nr,Pl,Pr,Bx,By);
    }


    std::vector<qk::solver::bc_periodic> electron_bcs;
    {
        qk::solver::bc_periodic bc;
        bc.add_output_variable(electrons[0]);
        electron_bcs.push_back(bc);
    }
    {
        qk::solver::bc_periodic bc;
        bc.add_output_variable(electrons[1]);
        electron_bcs.push_back(bc);
    }
    {
        qk::solver::bc_periodic bc;
        bc.add_output_variable(electrons[2]);
        electron_bcs.push_back(bc);
    }

    std::vector<qk::solver::bc_periodic> ion_bcs;
    {
        qk::solver::bc_periodic bc;
        bc.add_output_variable(ions[0]);
        ion_bcs.push_back(bc);
    }
    {
        qk::solver::bc_periodic bc;
        bc.add_output_variable(ions[1]);
        ion_bcs.push_back(bc);
    }
    {
        qk::solver::bc_periodic bc;
        bc.add_output_variable(ions[2]);
        ion_bcs.push_back(bc);
    }

    std::vector<qk::solver::bc_periodic> field_bcs;
    {
        qk::solver::bc_periodic bc;
        bc.add_output_variable(fields[0]);
        field_bcs.push_back(bc);
    }
    {
        qk::solver::bc_periodic bc;
        bc.add_output_variable(fields[1]);
        field_bcs.push_back(bc);
    }
    {
        qk::solver::bc_periodic bc;
        bc.add_output_variable(fields[2]);
        field_bcs.push_back(bc);
    }

    qk::solver::fill fill;
    {
        fill.setup(0.);
        fill.add_output_variable(electrons[4]);
        fill.add_output_variable(ions[4]);
        fill.add_output_variable(fields[4]);
    }

    qk::solver::swap swap;
    {
        swap.add_output_variable(electrons[0]);
        swap.add_output_variable(electrons[3]);
        swap.add_output_variable(ions[0]);
        swap.add_output_variable(ions[3]);
        swap.add_output_variable(fields[0]);
        swap.add_output_variable(fields[3]);
    }

    qk::variable::variable_manager variable_manager(grid);

    double average_delta_time = 0;

    ic.solve(variable_manager);

    std::vector<qk::variable::variable_id> write_vars;
    write_vars.push_back(electrons[0]);
    write_vars.push_back(ions[0]);
    write_vars.push_back(fields[0]);

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
            for(int k = 0; k < 3; k++){
                electron_bcs[k].solve(variable_manager);
                ion_bcs[k].solve(variable_manager);
                field_bcs[k].solve(variable_manager);
                fill.solve(variable_manager);
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
