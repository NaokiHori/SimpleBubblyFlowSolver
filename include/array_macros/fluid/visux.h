#if !defined(INCLUDE_ARRAY_MACROS_FLUID_VISUX_H)
#define INCLUDE_ARRAY_MACROS_FLUID_VISUX_H

// This file is generated by tools/define_arrays.py

// [1 : isize+1], [0 : jsize+1]
#define VISUX(I, J) (visux[(I-1) + (isize+1) * (J  )])
#define VISUX_NADDS (int [NDIMS][2]){ {0, 1}, {1, 1}, }

#endif // INCLUDE_ARRAY_MACROS_FLUID_VISUX_H
