#include <mpi.h>
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "interface.h"
#include "../internal.h"
#include "internal.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/normal.h"
#include "array_macros/interface/fluxx.h"

static int update_boundaries(
    const domain_t * domain,
    array_t * array
){
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    double * fluxx = array->data;
    // assume impermeable walls
    for(int j = 1; j <= jsize; j++){
      FLUXX(      1, j) = 0.;
      FLUXX(isize+1, j) = 0.;
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

int compute_flux_x(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict ux = fluid->ux.data;
  const double * restrict vof = interface->vof.data;
  const normal_t * restrict normal = interface->normal.data;
  double * restrict fluxx = interface->fluxx.data;
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      // use upwind information
      const double vel = UX(i, j);
      const int    ii = vel < 0. ?    i : i - 1;
      const double  x = vel < 0. ? -0.5 :  +0.5;
      // evaluate flux
      const double lvof = VOF(ii, j);
      if(lvof < vofmin || 1. - vofmin < lvof){
        FLUXX(i, j) = vel * lvof;
        continue;
      }
      double flux = 0.;
      for(int jj = 0; jj < NGAUSS; jj++){
        const double w = gauss_ws[jj];
        const double y = gauss_ps[jj];
        flux += w * indicator(NORMAL(ii, j), (const double [NDIMS]){x, y});
      }
      FLUXX(i, j) = vel * flux;
    }
  }
  update_boundaries(domain, &interface->fluxx);
  return 0;
}

