#include "ssprk3.h"

// STL include

// QK includes
#include "lib/exception.h"
#include "grid/rectilinear.h"
#include "variable/variable_id.h"
#include "variable/variable.h"

namespace qk
{
namespace solver
{

void ssprk3::solve(qk::variable::variable_manager & variable_manager, const int tag) const
{
    if (_input_variable_ids.size() != 2) {
        throw qk::exception("qk::solver::ssprk3::solve : Input must contain two variables (q_0, rhs).");
    }

    if (_output_variable_ids.size() < 3) {
        throw qk::exception("qk::solver::ssprk3::solve : Output must contain three variables (q_1, q_2, q_3).");
    }

    const double dt = variable_manager.dt();

    const qk::variable::variable_id & id_0 = _input_variable_ids[0];
    const qk::variable::variable & q_0 = variable_manager.input_variable(id_0);

    const qk::variable::variable_id & id_rhs = _input_variable_ids[1];
    const qk::variable::variable & rhs = variable_manager.input_variable(id_rhs);

    const qk::variable::variable_id & id_1 = _output_variable_ids[0];
    qk::variable::variable & q_1 = variable_manager.output_variable(id_1);

    const qk::variable::variable_id & id_2 = _output_variable_ids[1];
    qk::variable::variable & q_2 = variable_manager.output_variable(id_2);

    const qk::variable::variable_id & id_3 = _output_variable_ids[2];
    qk::variable::variable & q_3 = variable_manager.output_variable(id_3);

    for (qk::indexer chunk_idx = q_0.indexer(); chunk_idx.exists(); chunk_idx.next()) {
        if (tag == STAGE_0)
            stage_0(dt, q_0[chunk_idx], rhs[chunk_idx], q_1[chunk_idx]);

        if (tag == STAGE_1)
            stage_1(dt, q_0[chunk_idx], q_1[chunk_idx], rhs[chunk_idx], q_2[chunk_idx]);

        if (tag == STAGE_2)
            stage_2(dt, q_0[chunk_idx], q_1[chunk_idx], q_2[chunk_idx], rhs[chunk_idx], q_3[chunk_idx]);
    }

}

void ssprk3::stage_0(const double dt,
    const qk::data::extended_datachunk & q0,
    const qk::data::extended_datachunk & rhs,
    qk::data::extended_datachunk & q1) const
{
    const double * __restrict__ n0 = q0.data();
    const double * __restrict__ r = rhs.data();
    double * __restrict__ n1 = q1.data();
    const int num_points = rhs.volume();
    for (int i = 0; i < num_points; i++) {
        n1[i] = n0[i] + dt * r[i];
    }
}

void ssprk3::stage_1(const double dt,
    const qk::data::extended_datachunk & q0,
    const qk::data::extended_datachunk & q1,
    const qk::data::extended_datachunk & rhs,
    qk::data::extended_datachunk & q2) const
{
    const double * __restrict__ n0 = q0.data();
    const double * __restrict__ n1 = q1.data();
    const double * __restrict__ r = rhs.data();
    double * __restrict__ n2 = q2.data();
    const double c = 0.25 * dt;
    const int num_points = rhs.volume();
    for (int i = 0; i < num_points; i++) {
        n2[i] = 0.75 * n0[i] + 0.25 * n1[i] + c * r[i];
    }
}

void ssprk3::stage_2(const double dt,
    const qk::data::extended_datachunk & q0,
    const qk::data::extended_datachunk & q1,
    const qk::data::extended_datachunk & q2,
    const qk::data::extended_datachunk & rhs,
    qk::data::extended_datachunk & q3) const
{
    const double * __restrict__ n0 = q0.data();
    const double * __restrict__ n1 = q1.data();
    const double * __restrict__ n2 = q2.data();
    const double * __restrict__ r = rhs.data();
    double * __restrict__ n3 = q3.data();
    const int num_points = rhs.volume();
    for (int i = 0; i < num_points; i++) {
        n3[i] = (n0[i] + 2 * (n2[i] + dt * r[i])) / 3.;
    }
}

}
}
