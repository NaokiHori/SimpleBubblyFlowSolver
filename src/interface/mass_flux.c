#include "domain.h"
#include "fluid.h"
#include "interface.h"
#include "interface_solver.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/interface/fluxx.h"
#include "array_macros/interface/fluxy.h"

static int convert_x(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict ux = fluid->ux.data;
  double * restrict fluxx = interface->fluxx.data;
  const double denr = fluid->denr;
  // convert x vof flux to x mass flux
  for(int j = 0; j <= jsize + 1; j++){
    for(int i = 1; i <= isize + 1; i++){
      FLUXX(i, j) = UX(i, j) + (denr - 1.) * FLUXX(i, j);
    }
  }
  return 0;
}

static int convert_y(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict uy = fluid->uy.data;
  double * restrict fluxy = interface->fluxy.data;
  const double denr = fluid->denr;
  // convert y vof flux to y mass flux
  for(int j = 0; j <= jsize + 1; j++){
    for(int i = 0; i <= isize + 1; i++){
      FLUXY(i, j) = UY(i, j) + (denr - 1.) * FLUXY(i, j);
    }
  }
  return 0;
}

int interface_compute_mass_flux(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  // NOTE: since all data is inside the process,
  //         update halo and boundary values locally
  convert_x(domain, fluid, interface);
  convert_y(domain, fluid, interface);
  return 0;
}

