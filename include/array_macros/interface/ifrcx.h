#if !defined(INCLUDE_ARRAY_MACROS_INTERFACE_IFRCX_H)
#define INCLUDE_ARRAY_MACROS_INTERFACE_IFRCX_H

// This file is generated by tools/define_arrays.py

#if NDIMS == 2
// [2 : isize+0], [1 : jsize+0]
#define IFRCX(I, J) (ifrcx[(I-2) + (isize-1) * (J-1)])
#define IFRCX_NADDS (int [NDIMS][2]){ {-1, 0}, {0, 0}, }
#endif

#if NDIMS == 3
// [2 : isize+0], [1 : jsize+0], [1 : ksize+0]
#define IFRCX(I, J, K) (ifrcx[(I-2) + (isize-1) * ((J-1) + (jsize+0) * (K-1))])
#define IFRCX_NADDS (int [NDIMS][2]){ {-1, 0}, {0, 0}, {0, 0}, }
#endif

#endif // INCLUDE_ARRAY_MACROS_INTERFACE_IFRCX_H