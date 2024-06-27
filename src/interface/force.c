#include "domain.h"
#include "fluid.h"
#include "interface.h"
#include "internal.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/fluid/den.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/ifrcx.h"
#include "array_macros/interface/ifrcy.h"
#include "array_macros/interface/ifrcz.h"
#include "array_macros/interface/curv.h"

// compute density factor
static double compute_refdeninv(
    const double denr
){
  return 2. / (1. + denr);
}

static int compute_force_x(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double tension = 1. / interface->We;
  const double * restrict den = fluid->den[0].data;
  const double refdeninv = compute_refdeninv(fluid->denr);
  const double * restrict vof = interface->vof.data;
  const double * restrict curv = interface->curv.data;
  double * restrict ifrcx = interface->ifrcx.data;
  // compute surface tension force in x direction
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      const double den_x0 = (
          + 0.5 * DEN(i-1, j  )
          + 0.5 * DEN(i  , j  )
      );
      const double kappa = (
          + 0.5 * CURV(i-1, j  )
          + 0.5 * CURV(i  , j  )
      );
      const double delta = 1. / HXXF(i  ) * (
          - VOF(i-1, j  )
          + VOF(i  , j  )
      );
      IFRCX(i, j) = den_x0 * refdeninv * tension * kappa * delta;
    }
  }
  return 0;
}

static int compute_force_y(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double hy = domain->hy;
  const double tension = 1. / interface->We;
  const double * restrict den = fluid->den[0].data;
  const double refdeninv = compute_refdeninv(fluid->denr);
  const double * restrict vof = interface->vof.data;
  const double * restrict curv = interface->curv.data;
  double * restrict ifrcy = interface->ifrcy.data;
  // compute surface tension force in y direction
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      const double den_y0 = (
          + 0.5 * DEN(i  , j-1)
          + 0.5 * DEN(i  , j  )
      );
      const double kappa = (
          + 0.5 * CURV(i  , j-1)
          + 0.5 * CURV(i  , j  )
      );
      const double delta = 1. / hy * (
          - VOF(i  , j-1)
          + VOF(i  , j  )
      );
      IFRCY(i, j) = den_y0 * refdeninv * tension * kappa * delta;
    }
  }
  return 0;
}

int interface_compute_force(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  compute_force_x(domain, fluid, interface);
  compute_force_y(domain, fluid, interface);
  return 0;
}

