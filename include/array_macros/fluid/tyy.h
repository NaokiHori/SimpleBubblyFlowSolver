#if !defined(INCLUDE_ARRAY_MACROS_FLUID_TYY_H)
#define INCLUDE_ARRAY_MACROS_FLUID_TYY_H

// This file is generated by tools/define_arrays.py

#if NDIMS == 2
// [1 : isize+0], [0 : jsize+1]
#define TYY(I, J) (tyy[(I-1) + (isize+0) * (J  )])
#define TYY_NADDS (int [NDIMS][2]){ {0, 0}, {1, 1}, }
#endif

#if NDIMS == 3
// [1 : isize+0], [0 : jsize+1], [0 : ksize+1]
#define TYY(I, J, K) (tyy[(I-1) + (isize+0) * ((J  ) + (jsize+2) * (K  ))])
#define TYY_NADDS (int [NDIMS][2]){ {0, 0}, {1, 1}, {1, 1}, }
#endif

#endif // INCLUDE_ARRAY_MACROS_FLUID_TYY_H
