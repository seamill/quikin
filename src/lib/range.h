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

    bool includes(const qk::range & rng) const;

    void set(const int dim, const int lower, const int upper);

    range & operator=(const range & range);

    friend std::ostream &
    operator<<(std::ostream & os, const qk::range & range)
    {
        os << "Range : lower [";
        for(int i = 0; i < range.num_dims(); ++i){
            os << range.lower(i);
            if(i != range.num_dims()-1){
                os << ", ";
            }
        }
        os <<"], upper [";
        for(int i = 0; i < range.num_dims(); ++i){
            os << range.upper(i);
            if(i != range.num_dims()-1){
                os << ", ";
            }
        }
        os<<"]";
    }

protected:

    // This is the only function with write access to _num_dims, _lower, _upper_, and _stride.
    void setup(const int numDims, const int *lower, const int *upper);

    int _num_dims;
    std::vector<int> _lower;
    std::vector<int> _upper;
    std::vector<int> _stride;

private:

};

}


inline bool operator==(const qk::range & a, const qk::range & b)
{

    if(a.num_dims()!=b.num_dims()){
        return false;
    }

    for(int i=0;i<a.num_dims();++i){
        if(a.lower(i) != b.lower(i)){
            return false;
        }
        if(a.upper(i) != b.upper(i)){
            return false;
        }
    }

    return true;
}

inline bool operator!=(const qk::range & a, const qk::range & b)
{
    return !(a==b);
}
#endif // _qk_lib_range_H
