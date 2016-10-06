
#include "test_interface.h"

// STL includes
#include <vector>

// QK includes
#include "variable/variable.h"
#include "basis/basis.h"
#include "basis/volume_average.h"

bool
test_variable(qk::test::test_monitor & tm)
{
    {
        qk::variable::variable var;
        tm.assert(var.volume() == 0, "qk::test::test_variable : Empty Constructor.");
    }

    {

        const std::string var_name = "CoMpLeX VaRiAbLe NaMe";

        std::vector<std::string> comp_names;
        comp_names.push_back("A");
        comp_names.push_back("B");
        comp_names.push_back("C");
        qk::basis::basis basis = qk::basis::volume_average();
        int num_dims = 4;
        int dims[7] = {7,4,3,1};
        const qk::range rng(4,dims);
        const int num_ghost_layers = 3;

        qk::range ext_rng(rng);
        for(int i = 0; i < num_dims; ++i){
            ext_rng.expand(i,-num_ghost_layers,num_ghost_layers);
        }

        qk::variable::variable var(var_name,comp_names, basis, rng,num_ghost_layers);

        tm.assert(var.name() == var_name, "qk::test::test_variable : Main constructor name test.");
        tm.assert(var.component_names() == comp_names, "qk::test::test_variable : Main constructor component names test.");

        int num = 0;
        for(qk::indexer idx = var.indexer(); idx.exists(); ++idx){
            const qk::data::extended_datachunk & dchunk = var[idx];

            tm.assert(dchunk.volume() == ext_rng.volume()*comp_names.size()*basis.num_points(), "qk::test::test_variable : Main constructor volume test.");
            tm.assert(dchunk.length(num_dims) == basis.num_points(), "qk::test::test_variable : Main constructor basis range test.");
            tm.assert(dchunk.length(num_dims+1) == comp_names.size(), "qk::test::test_variable : Main constructor component range test.");

            ++num;
        }
        tm.assert(num==1, "qk::test::test_variable : Main constructor number of datachunks test.");

    }

    return true;
}

QK_TEST_CALL_GROUP_START("Variable tests")
QK_TEST_CALL("Some variable tests", test_variable);
QK_TEST_CALL_GROUP_END("Variable tests")
