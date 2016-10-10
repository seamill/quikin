#include <omp.h>

#define _QK_NUM_THREADS_ 8

#define _QK_DOUBLE_ALIGNMENT_ 64
#define _QK_ASSUME_ALIGNED(ptr) (double*)__builtin_assume_aligned(ptr,_QK_DOUBLE_ALIGNMENT_)
#define _QK_ASSUME_ALIGNED_CONST(ptr) (const double*)__builtin_assume_aligned(ptr,_QK_DOUBLE_ALIGNMENT_)

