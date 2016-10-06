
#include "test_interface.h"

// QK includes
#include "lib/indexer.h"
#include "lib/indexer_interface.h"


#include <iostream>
#include <iomanip>

bool
test_constructors(qk::test::test_monitor & tm)
{
    {
        qk::indexer_interface<int> ii;
        tm.assert(ii.volume() == 0, "qk::test::test_indexer : Empty Indexer Interface Constructor.");
    }

    {
        int num_dims = 7;
        int dims[7] = {7,4,3,1,8,9,2};
        qk::range rng(7,dims);

        const qk::indexer_interface<int> ii(rng);

        tm.assert(ii.num_dims()==7, "qk::test::test_indexer : Indexer Interface range construction test.");

        const qk::indexer_interface<int> ii2(ii);

        tm.assert(ii2.num_dims()==7, "qk::test::test_indexer : Indexer Interface copy construction test.");

    }

    return true;
}

bool
test_functionality(qk::test::test_monitor & tm)
{

    qk::range rng;
    rng.extrude(-2,4);
    rng.extrude(-5,90);

    {
        qk::range rng1;
        rng1.extrude(0,100);

        qk::indexer_interface<int> ii(rng);

        ii.resize(rng1);

        tm.assert(ii.volume() == 100, "qk::test::test_indexer : Indexer Interface resize test.");

    }

    {
        qk::indexer_interface<int> ii(rng);

        ii.fill(17);

        tm.assert((ii.data()[0] == 17) && (ii.data()[ii.volume()-1] == 17), "qk::test::test_indexer : Fill test.");

        qk::indexer_interface<int> ii2(rng);

        ii2.swap(ii);

        tm.assert((ii2.data()[0] == 17) && (ii2.data()[ii2.volume()-1] == 17), "qk::test::test_indexer : Swap test.");

    }

    return true;
}


bool
test_indexer(qk::test::test_monitor & tm)
{

    // Create a 100 by 100 array
    qk::range rng;
    rng.extrude(0,15);
    rng.extrude(0,25);
    qk::indexer_interface<int> ii(rng);

    qk::range sub_rng;
    sub_rng.extrude(9,13);
    sub_rng.extrude(3,7);

    qk::range fill_rng;
    fill_rng.extrude(1,5);
    fill_rng.extrude(9,13);

    {
        qk::indexer idx;
        tm.assert(! idx.exists(), "qk::test::test_indexer : Null Indexer not exists test.");
    }

    {
        qk::indexer idx = qk::indexer(sub_rng, rng);
        tm.assert(idx.exists(), "qk::test::test_indexer : Sub-range Indexer exists test.");

        tm.assert(idx[0] == sub_rng.lower(0), "qk::test::test_indexer : Sub-range Indexer dim index test 0.");
        tm.assert(idx[1] == sub_rng.lower(1), "qk::test::test_indexer : Sub-range Indexer dim index test 1.");

    }

    {
        qk::indexer idx = ii.indexer(qk::range());
        tm.assert(! idx.exists(), "qk::test::test_indexer : Null Indexer from Indexer Interface exists test.");
    }

    {
        qk::indexer idx = ii.indexer();
        tm.assert(idx.exists(), "qk::test::test_indexer : Indexer from Indexer Interface exists test.");
        tm.assert(idx.linear_index()==0, "qk::test::test_indexer : Indexer linear index initialized to zero test.");
        tm.assert(idx[0]==ii.lower(0), "qk::test::test_indexer : Indexer dim index initialized to lower test 1.");
        tm.assert(idx[1]==ii.lower(1), "qk::test::test_indexer : Indexer dim index initialized to lower test 2.");
    }

    {
        int num = 0;
        for(qk::indexer idx = ii.indexer();idx.exists();idx.next()){
            ++num;
        }
        tm.assert(num==rng.volume(), "qk::test::test_indexer : Indexer full iteration test.");
    }

    {
        int num = 0;
        for(qk::indexer idx = ii.indexer(sub_rng);idx.exists();idx.next()){
            ++num;
        }
        tm.assert(num==sub_rng.volume(), "qk::test::test_indexer : Indexer partial iteration test.");
    }

    {
        qk::range rng3;
        rng.extrude(-1,39);
        qk::indexer idx = ii.indexer(rng3);
        tm.assert(! idx.exists(), "qk::test::test_indexer : Dimensional mismatch Indexer request from Indexer Interface does not exist test.");
    }

    {
        qk::range rng3;
        rng.extrude(-100,100);
        rng.extrude(-100,100);
        qk::indexer idx = ii.indexer(rng3);
        tm.assert(!idx.exists(), "qk::test::test_indexer : Super-range Indexer request from Indexer Interface does not exist test.");
    }

    {


        {
            qk::indexer idx = ii.indexer(fill_rng);
            tm.assert(idx[0]==fill_rng.lower(0), "qk::test::test_indexer : Indexer request from Indexer Interface dim index initialized to lower of request test 1.");
            tm.assert(idx[1]==fill_rng.lower(1), "qk::test::test_indexer : Indexer request from Indexer Interface dim index initialized to lower of request test 2.");
        }


        ii.fill(0);

        {

            for(qk::indexer idx = ii.indexer(fill_rng); idx.exists();++idx){
                ii[idx] = 21;
            }

            bool test_21_value = true;
            bool test_0_value = true;
            for(qk::indexer idx = ii.indexer(); idx.exists();++idx){
                const bool in_21_range = ((idx[0] >= fill_rng.lower(0)) and (idx[0] < fill_rng.upper(0))) and ((idx[1] >= fill_rng.lower(1)) and (idx[1] < fill_rng.upper(1)));

                if(in_21_range){
                    test_21_value = test_21_value && (ii[idx] == 21);
                } else {
                    test_0_value = test_0_value && (ii[idx] == 0);
                }
            }
            tm.assert(test_0_value, "qk::test::test_indexer : Sub-range coloring red test.");
            tm.assert(test_21_value, "qk::test::test_indexer : Sub-range coloring black test.");

        }

        qk::indexer_interface<int> ii2(sub_rng);


        {
            ii2.pull(ii,fill_rng,ii2.range());
            bool test_21_value = true;
            for(qk::indexer idx = ii2.indexer(); idx.exists();++idx){
                test_21_value = test_21_value && (ii2[idx] == 21);
            }
            tm.assert(test_21_value, "qk::test::test_indexer : Sub-range pull test.");
        }

        {
            ii2.fill(7);
            ii.pull(ii2,ii2.range(),fill_rng);

            bool test_7_value = true;
            bool test_0_value = true;
            for(qk::indexer idx = ii.indexer(); idx.exists();++idx){
                const bool in_7_range = ((idx[0] >= fill_rng.lower(0)) and (idx[0] < fill_rng.upper(0))) and ((idx[1] >= fill_rng.lower(1)) and (idx[1] < fill_rng.upper(1)));

                if(in_7_range){
                    test_7_value = test_7_value && (ii[idx] == 7);
                } else {
                    test_0_value = test_0_value && (ii[idx] == 0);
                }
            }

            tm.assert(test_0_value, "qk::test::test_indexer : Sub-range pull coloring red test.");
            tm.assert(test_7_value, "qk::test::test_indexer : Sub-range pull coloring black test.");

        }

        {
            ii.fill(fill_rng, 15);

            bool test_15_value = true;
            bool test_0_value = true;
            for(qk::indexer idx = ii.indexer(); idx.exists();++idx){
                const bool in_15_range = ((idx[0] >= fill_rng.lower(0)) and (idx[0] < fill_rng.upper(0))) and ((idx[1] >= fill_rng.lower(1)) and (idx[1] < fill_rng.upper(1)));

                if(in_15_range){
                    test_15_value = test_15_value && (ii[idx] == 15);
                } else {
                    test_0_value = test_0_value && (ii[idx] == 0);
                }
            }

            tm.assert(test_0_value, "qk::test::test_indexer : Sub-range fill coloring red test.");
            tm.assert(test_15_value, "qk::test::test_indexer : Sub-range fill coloring black test.");

        }

    }

    return true;
}

QK_TEST_CALL_GROUP_START("Indexer tests")
QK_TEST_CALL("Indexer Interface constructor tests", test_constructors);
QK_TEST_CALL("Indexer Interface functionality tests", test_functionality);
QK_TEST_CALL("Indexer Interface indexer tests", test_indexer);
QK_TEST_CALL_GROUP_END("Indexer tests")
