#if !defined(INTERFACE_H)
#define INTERFACE_H

#include "array.h"

// data type for N-dimensional vectors
typedef double vector_t[NDIMS];
// data type for (N+1)-dimensional vectors
//   to store additional element (segment),
//   which is stored at the last
typedef double normal_t[NDIMS + 1];

typedef struct {
  array_t vof;
  array_t ifrcx;
  array_t ifrcy;
#if NDIMS == 3
  array_t ifrcz;
#endif
  array_t dvof;
  array_t normal;
  array_t curv;
  array_t fluxx;
  array_t fluxy;
#if NDIMS == 3
  array_t fluxz;
#endif
  array_t src[2];
  double We;
} interface_t;

#endif // INTERFACE_H
