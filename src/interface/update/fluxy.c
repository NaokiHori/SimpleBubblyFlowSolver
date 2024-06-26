#include <mpi.h>
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "interface.h"
#include "../internal.h"
#include "internal.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/normal.h"
#include "array_macros/interface/fluxy.h"

static int update_boundaries(
    const domain_t * domain,
    array_t * array
){
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    double * fluxy = array->data;
    // dummy
    for(int j = 1; j <= jsize; j++){
      FLUXY(      0, j) = 0.;
      FLUXY(isize+1, j) = 0.;
    }
  }
  {
    static MPI_Datatype dtypes[NDIMS - 1] = {
      MPI_DOUBLE,
    };
    if(0 != halo_communicate_in_y(domain, dtypes + 0, array)){
      return 1;
    }
  }
  return 0;
}

int compute_flux_y(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict uy = fluid->uy.data;
  const double * restrict vof = interface->vof.data;
  const normal_t * restrict normal = interface->normal.data;
  double * restrict fluxy = interface->fluxy.data;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      // use upwind information
      const double vel = UY(i, j);
      const int    jj = vel < 0. ?    j : j - 1;
      const double  y = vel < 0. ? -0.5 :  +0.5;
      // evaluate flux
      const double lvof = VOF(i, jj);
      if(lvof < vofmin || 1. - vofmin < lvof){
        FLUXY(i, j) = vel * lvof;
        continue;
      }
      double flux = 0.;
      for(int ii = 0; ii < NGAUSS; ii++){
        const double w = gauss_ws[ii];
        const double x = gauss_ps[ii];
        flux += w * indicator(NORMAL(i, jj), (const double [NDIMS]){x, y});
      }
      FLUXY(i, j) = vel * flux;
    }
  }
  update_boundaries(domain, &interface->fluxy);
  return 0;
}

