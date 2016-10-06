
#include "test_interface.h"

// QK includes
#include "lib/range.h"

bool
test_constructors(qk::test::test_monitor & tm)
{
    {
        qk::range rng;
        tm.assert(rng.volume() == 0, "qk::test::test_range : Empty Constructor.");
    }

    {
        int num_dims = 7;
        int dims[7] = {7,4,3,1,8,9,2};
        qk::range rng(7,dims);

        int stride[num_dims];
        stride[num_dims-1]=1;
        for(int i = num_dims-1; i > 0; --i){
            stride[i-1] = stride[i] * dims[i];
        }

        tm.assert(rng.num_dims()==7, "qk::test::test_range : Block constructor num_dims test.");

        bool lower_test = true, upper_test = true, length_test = true, stride_test = true;
        for(int i = 0; i < 7; ++i){
            lower_test = lower_test && (rng.lower(i) == 0);
            upper_test = upper_test && (rng.upper(i) == dims[i]);
            length_test = length_test && (rng.length(i) == dims[i]);
            stride_test = stride_test && (rng.stride(i) == stride[i]);
        }
        tm.assert(lower_test, "qk::test::test_range : Block constructor lower test.");
        tm.assert(upper_test, "qk::test::test_range : Block constructor upper test.");
        tm.assert(length_test, "qk::test::test_range : Block constructor length test.");
        tm.assert(stride_test, "qk::test::test_range : Block constructor stride test.");

        tm.assert(rng.volume() == 7*4*3*1*8*9*2, "qk::test::test_range : Block constructor volume test.");

    }

    {
        int num_dims=3;
        int lower[3] = {1,-20,8};
        int upper[3] = {7,-10,9};
        int stride[num_dims];
        stride[num_dims-1]=1;
        for(int i = num_dims-1; i > 0; --i){
            stride[i-1] = stride[i] * (upper[i] - lower[i]);
        }


        qk::range rng(3,lower,upper);

        tm.assert(rng.num_dims()==3, "qk::test::test_range : Bound constructor num_dims test.");

        {
            bool lower_test = true, upper_test = true, length_test = true, stride_test = true;
            for(int i = 0; i < 3; ++i){
                lower_test = lower_test && (rng.lower(i) == lower[i]);
                upper_test = upper_test && (rng.upper(i) == upper[i]);
                length_test = length_test && (rng.length(i) == upper[i] - lower[i]);
                stride_test = stride_test && (rng.stride(i) == stride[i]);
            }
            tm.assert(lower_test, "qk::test::test_range : Bound constructor lower test.");
            tm.assert(upper_test, "qk::test::test_range : Bound constructor upper test.");
            tm.assert(length_test, "qk::test::test_range : Bound constructor length test.");
            tm.assert(stride_test, "qk::test::test_range : Bound constructor stride test.");
        }

        qk::range rng2(rng);

        {
            bool lower_test = true, upper_test = true, length_test = true, stride_test = true;
            for(int i = 0; i < 3; ++i){
                lower_test = lower_test && (rng2.lower(i) == lower[i]);
                upper_test = upper_test && (rng2.upper(i) == upper[i]);
                length_test = length_test && (rng2.length(i) == upper[i] - lower[i]);
                stride_test = stride_test && (rng2.stride(i) == stride[i]);
            }
            tm.assert(lower_test, "qk::test::test_range : Copy constructor lower test.");
            tm.assert(upper_test, "qk::test::test_range : Copy constructor upper test.");
            tm.assert(length_test, "qk::test::test_range : Copy constructor length test.");
            tm.assert(stride_test, "qk::test::test_range : Copy constructor stride test.");
        }

    }

    return true;
}

bool
test_extrusion(qk::test::test_monitor & tm)
{

    int num_dims = 3;
    int dims[3] = {10,10,10};

    const qk::range rng_base(num_dims,dims);

    // Test extrusion setup process
    {
        qk::range rng;

        rng.extrude(0,10);
        rng.extrude(0,10);
        rng.extrude(0,10);

        tm.assert(rng == rng_base, "qk::test::test_range : Extrude empty range.");
    }

    // Test extrusion process
    {
        qk::range rng(rng_base);

        rng.extrude(-51,-41);

        tm.assert(rng.volume() == rng_base.volume()*10, "qk::test::test_range : Extrude existing range.");
    }

    return true;
}

bool
test_resize(qk::test::test_monitor & tm)
{
    int num_dims = 3;
    int dims[3] = {10,10,10};

    const qk::range rng_base(num_dims,dims);

    // Test resize of existing range by setting lower and upper
    {
        qk::range rng(rng_base);

        int lower[2] = {2,7};
        int upper[2] = {5,9};
        rng.resize(2,lower,upper);

        tm.assert(rng.num_dims()==2, "qk::test::test_range : Resize num_dims test.");

        bool lower_test = true, upper_test = true;
        for(int i = 0; i < 2; ++i){
            lower_test = lower_test && (rng.lower(i) == lower[i]);
            upper_test = upper_test && (rng.upper(i) == upper[i]);

        }
        tm.assert(lower_test, "qk::test::test_range : Resize lower test.");
        tm.assert(upper_test, "qk::test::test_range : Resize upper test.");

        tm.assert(rng.volume() == 6, "qk::test::test_range : Resize volume test.");
    }

    // Test resize of existing range with another range
    {
        int lower[2] = {-2,90};
        int upper[2] = {48,982};
        const qk::range new_rng(2,lower, upper);

        qk::range rng(rng_base);

        rng.resize(new_rng);

        tm.assert(rng == new_rng, "qk::test::test_range : Resize from existing range.");
    }

    return true;
}

bool
test_functionality(qk::test::test_monitor & tm)
{
    int num_dims = 3;
    int dims[3] = {10,10,10};

    const qk::range rng_base(num_dims,dims);

    {
        qk::range rng(rng_base);

        rng.shift(1,5);

        tm.assert(rng != rng_base, "qk::test::test_range : Shift did something test.");
        tm.assert(rng.lower(1) == 5, "qk::test::test_range : Shift lower test.");
        tm.assert(rng.upper(1) == 15, "qk::test::test_range : Shift upper test.");
        tm.assert(rng.stride(0) == rng_base.stride(0), "qk::test::test_range : Shift stride test.");
        tm.assert(rng.volume() == rng_base.volume(), "qk::test::test_range : Shift volume test.");

    }

    {
        qk::range rng(rng_base);

        rng.set(1,5,15);

        tm.assert(rng != rng_base, "qk::test::test_range : Set did something test.");
        tm.assert(rng.lower(1) == 5, "qk::test::test_range : Set lower test.");
        tm.assert(rng.upper(1) == 15, "qk::test::test_range : Set upper test.");
        tm.assert(rng.stride(0) == rng_base.stride(0), "qk::test::test_range : Set stride test.");
        tm.assert(rng.volume() == rng_base.volume(), "qk::test::test_range : Set volume test.");

    }

    {
        qk::range rng(rng_base);

        qk::range in_rng;
        in_rng.extrude(2,4);
        in_rng.extrude(2,4);

        tm.assert(!rng.includes(in_rng), "qk::test::test_range : Include sub-range dimension mismatch fail test.");

        in_rng.extrude(2,4);
        tm.assert(rng.includes(in_rng), "qk::test::test_range : Include sub-range test.");

        in_rng.set(1,-1, 11);
        tm.assert(!rng.includes(in_rng), "qk::test::test_range : Include sub-range oversized fail test.");

    }

    return true;

}

QK_TEST_CALL_GROUP_START("Range tests")
QK_TEST_CALL("Range constructor tests", test_constructors);
QK_TEST_CALL("Range extrusion tests", test_extrusion);
QK_TEST_CALL("Range resize tests", test_resize);
QK_TEST_CALL("Range functionality tests", test_functionality);
QK_TEST_CALL_GROUP_END("Range tests")
