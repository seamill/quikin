#ifndef _qk_lib_range_H
#define _qk_lib_range_H

// STL includes
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <cstring>

// This is a non inclusive range
namespace qk
{

class range
{
public:
    range();
    range(const range & range);
    range(const int numDims, const int * dims);
    range(const int numDims, const int * lower, const int * upper);

    virtual ~range();

    virtual void resize(const range & range);
    virtual void resize(const int numDims, const int *lower, const int *upper);

    void extrude(const int lower_extrude, const int upper_extrude);
    void expand(const int dim, const int lowerExpand, const int upperExpand);
    void shift(const int dim, const int shift);

    int num_dims() const;
    int length(const int i) const;
    int stride(const int i) const;

    int volume() const;

    int lower(const int dim) const {return _lower[dim];}
    int upper(const int dim) const {return _upper[dim];}

    void set(const int dim, const int lower, const int upper);

    range & operator=(const range & range);

protected:

    void setup(const int numDims, const int *lower, const int *upper);

    int _num_dims;
    std::vector<int> _lower;
    std::vector<int> _upper;
    std::vector<int> _stride;

private:

};

}

#endif // _qk_lib_range_H
