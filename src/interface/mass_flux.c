#include "domain.h"
#include "fluid.h"
#include "interface.h"
#include "interface_solver.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#if NDIMS == 3
#include "array_macros/fluid/uz.h"
#endif
#include "array_macros/interface/fluxx.h"
#include "array_macros/interface/fluxy.h"
#if NDIMS == 3
#include "array_macros/interface/fluxz.h"
#endif

static int convert_x(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict ux = fluid->ux.data;
  double * restrict fluxx = interface->fluxx.data;
  const double denr = fluid->denr;
#if NDIMS == 2
  // convert x vof flux to x mass flux | 5
  for(int j = 0; j <= jsize + 1; j++){
    for(int i = 1; i <= isize + 1; i++){
      FLUXX(i, j) = UX(i, j) + (denr - 1.) * FLUXX(i, j);
    }
  }
#else
  // convert x vof flux to x mass flux | 7
  for(int k = 0; k <= ksize + 1; k++){
    for(int j = 0; j <= jsize + 1; j++){
      for(int i = 1; i <= isize + 1; i++){
        FLUXX(i, j, k) = UX(i, j, k) + (denr - 1.) * FLUXX(i, j, k);
      }
    }
  }
#endif
  return 0;
}

static int convert_y(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict uy = fluid->uy.data;
  double * restrict fluxy = interface->fluxy.data;
  const double denr = fluid->denr;
#if NDIMS == 2
  // convert y vof flux to y mass flux | 5
  for(int j = 0; j <= jsize + 1; j++){
    for(int i = 0; i <= isize + 1; i++){
      FLUXY(i, j) = UY(i, j) + (denr - 1.) * FLUXY(i, j);
    }
  }
#else
  // convert y vof flux to y mass flux | 7
  for(int k = 0; k <= ksize + 1; k++){
    for(int j = 0; j <= jsize + 1; j++){
      for(int i = 0; i <= isize + 1; i++){
        FLUXY(i, j, k) = UY(i, j, k) + (denr - 1.) * FLUXY(i, j, k);
      }
    }
  }
#endif
  return 0;
}

#if NDIMS == 3
static int convert_z(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict uz = fluid->uz.data;
  double * restrict fluxz = interface->fluxz.data;
  const double denr = fluid->denr;
  // convert z vof flux to z mass flux | 7
  for(int k = 0; k <= ksize + 1; k++){
    for(int j = 0; j <= jsize + 1; j++){
      for(int i = 0; i <= isize + 1; i++){
        FLUXZ(i, j, k) = UZ(i, j, k) + (denr - 1.) * FLUXZ(i, j, k);
      }
    }
  }
  return 0;
}
#endif

int interface_compute_mass_flux(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  // NOTE: since all data is inside the process,
  //         update halo and boundary values locally
  convert_x(domain, fluid, interface);
  convert_y(domain, fluid, interface);
#if NDIMS == 3
  convert_z(domain, fluid, interface);
#endif
  return 0;
}

