#if !defined(INCLUDE_ARRAY_MACROS_INTERFACE_FLUXZ_H)
#define INCLUDE_ARRAY_MACROS_INTERFACE_FLUXZ_H

// This file is generated by tools/define_arrays.py

// [0 : isize+1], [0 : jsize+1], [0 : ksize+1]
#define FLUXZ(I, J, K) (fluxz[(I  ) + (isize+2) * ((J  ) + (jsize+2) * (K  ))])
#define FLUXZ_NADDS (int [NDIMS][2]){ {1, 1}, {1, 1}, {1, 1}, }

#endif // INCLUDE_ARRAY_MACROS_INTERFACE_FLUXZ_H
