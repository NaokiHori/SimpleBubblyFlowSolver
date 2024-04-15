#include "internal.h"

// pm 1 / 2 / sqrt(3)
const double gauss_ps[NGAUSS] = {
  - 0.2886751345948129,
  + 0.2886751345948129,
};
const double gauss_ws[NGAUSS] = {
  + 0.5,
  + 0.5,
};

const double vofbeta = 1.;
const double vofmin = 1.e-8;

