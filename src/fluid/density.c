#include <math.h>
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "interface.h"
#include "array_macros/fluid/den.h"
#include "array_macros/interface/vof.h"

int fluid_compute_density(
    const domain_t * domain,
    fluid_t * fluid,
    const interface_t * interface,
    const size_t index
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  double * restrict den = fluid->den[index].data;
  const double * restrict vof = interface->vof.data;
  const double denr = fluid->denr;
  const double min = fmin(1., denr);
  const double max = fmax(1., denr);
  for(int j = 0; j <= jsize + 1; j++){
    for(int i = 0; i <= isize + 1; i++){
      double * lden = &DEN(i, j);
      *lden = 1. + (denr - 1.) * VOF(i, j);
      *lden = fmax(min, *lden);
      *lden = fmin(max, *lden);
    }
  }
  return 0;
}

