#if !defined(INCLUDE_ARRAY_MACROS_STATISTICS_VOF1_H)
#define INCLUDE_ARRAY_MACROS_STATISTICS_VOF1_H

// This file is generated by tools/define_arrays.py

// [0 : isize+1], [1 : jsize+0]
#define VOF1(I, J) (vof1[(I  ) + (isize+2) * (J-1)])
#define VOF1_NADDS (int [NDIMS][2]){ {1, 1}, {0, 0}, }

#endif // INCLUDE_ARRAY_MACROS_STATISTICS_VOF1_H
