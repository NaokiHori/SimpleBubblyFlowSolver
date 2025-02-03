#if !defined(INCLUDE_ARRAY_MACROS_FLUID_PSI_H)
#define INCLUDE_ARRAY_MACROS_FLUID_PSI_H

// This file is generated by tools/define_arrays.py

// [0 : isize+1], [0 : jsize+1]
#define PSI(I, J) (psi[(I  ) + (isize+2) * (J  )])
#define PSI_NADDS (int [NDIMS][2]){ {1, 1}, {1, 1}, }

#endif // INCLUDE_ARRAY_MACROS_FLUID_PSI_H
