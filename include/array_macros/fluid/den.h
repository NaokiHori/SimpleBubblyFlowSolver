#if !defined(INCLUDE_ARRAY_MACROS_FLUID_DEN_H)
#define INCLUDE_ARRAY_MACROS_FLUID_DEN_H

// This file is generated by tools/define_arrays.py

#if NDIMS == 2
// [0 : isize+1], [0 : jsize+1]
#define DEN(I, J) (den[(I  ) + (isize+2) * (J  )])
#define DEN_NADDS (int [NDIMS][2]){ {1, 1}, {1, 1}, }
#endif

#if NDIMS == 3
// [0 : isize+1], [0 : jsize+1], [0 : ksize+1]
#define DEN(I, J, K) (den[(I  ) + (isize+2) * ((J  ) + (jsize+2) * (K  ))])
#define DEN_NADDS (int [NDIMS][2]){ {1, 1}, {1, 1}, {1, 1}, }
#endif

#endif // INCLUDE_ARRAY_MACROS_FLUID_DEN_H