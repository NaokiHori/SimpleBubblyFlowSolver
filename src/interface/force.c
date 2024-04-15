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

// compute density factor | 5
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
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
  const double tension = 1. / interface->We;
  const double * restrict den = fluid->den[0].data;
  const double refdeninv = compute_refdeninv(fluid->denr);
  const double * restrict vof = interface->vof.data;
  const double * restrict curv = interface->curv.data;
  double * restrict ifrcx = interface->ifrcx.data;
#if NDIMS == 2
  // compute surface tension force in x direction | 17
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
#else
  // compute surface tension force in x direction | 19
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        const double den_x0 = (
            + 0.5 * DEN(i-1, j  , k  )
            + 0.5 * DEN(i  , j  , k  )
        );
        const double kappa = (
            + 0.5 * CURV(i-1, j  , k  )
            + 0.5 * CURV(i  , j  , k  )
        );
        const double delta = 1. / HXXF(i  ) * (
            - VOF(i-1, j  , k  )
            + VOF(i  , j  , k  )
        );
        IFRCX(i, j, k) = den_x0 * refdeninv * tension * kappa * delta;
      }
    }
  }
#endif
  return 0;
}

static int compute_force_y(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double hy = domain->hy;
  const double tension = 1. / interface->We;
  const double * restrict den = fluid->den[0].data;
  const double refdeninv = compute_refdeninv(fluid->denr);
  const double * restrict vof = interface->vof.data;
  const double * restrict curv = interface->curv.data;
  double * restrict ifrcy = interface->ifrcy.data;
#if NDIMS == 2
  // compute surface tension force in y direction | 17
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
#else
  // compute surface tension force in y direction | 19
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double den_y0 = (
            + 0.5 * DEN(i  , j-1, k  )
            + 0.5 * DEN(i  , j  , k  )
        );
        const double kappa = (
            + 0.5 * CURV(i  , j-1, k  )
            + 0.5 * CURV(i  , j  , k  )
        );
        const double delta = 1. / hy * (
            - VOF(i  , j-1, k  )
            + VOF(i  , j  , k  )
        );
        IFRCY(i, j, k) = den_y0 * refdeninv * tension * kappa * delta;
      }
    }
  }
#endif
  return 0;
}

#if NDIMS == 3
static int compute_force_z(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double tension = 1. / interface->We;
  const double * restrict den = fluid->den[0].data;
  const double refdeninv = compute_refdeninv(fluid->denr);
  const double * restrict vof = interface->vof.data;
  const double * restrict curv = interface->curv.data;
  double * restrict ifrcz = interface->ifrcz.data;
  // compute surface tension force in z direction | 19
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double den_z0 = (
            + 0.5 * DEN(i  , j  , k-1)
            + 0.5 * DEN(i  , j  , k  )
        );
        const double kappa = (
            + 0.5 * CURV(i  , j  , k-1)
            + 0.5 * CURV(i  , j  , k  )
        );
        const double delta = 1. / hz * (
            - VOF(i  , j  , k-1)
            + VOF(i  , j  , k  )
        );
        IFRCZ(i, j, k) = den_z0 * refdeninv * tension * kappa * delta;
      }
    }
  }
  return 0;
}
#endif

int interface_compute_force(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  compute_force_x(domain, fluid, interface);
  compute_force_y(domain, fluid, interface);
#if NDIMS == 3
  compute_force_z(domain, fluid, interface);
#endif
  return 0;
}

