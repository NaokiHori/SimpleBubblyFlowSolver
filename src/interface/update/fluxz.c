#include <mpi.h>
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "interface.h"
#include "../internal.h"
#include "internal.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/normal.h"
#include "array_macros/interface/fluxz.h"

static int update_boundaries(
    const domain_t * domain,
    array_t * array
){
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const int ksize = domain->mysizes[2];
    double * fluxz = array->data;
    // dummy
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        FLUXZ(      0, j, k) = 0.;
        FLUXZ(isize+1, j, k) = 0.;
      }
    }
  }
  {
    static MPI_Datatype dtypes[NDIMS - 1] = {
      MPI_DOUBLE,
      MPI_DOUBLE,
    };
    if(0 != halo_communicate_in_y(domain, dtypes + 0, array)){
      return 1;
    }
    if(0 != halo_communicate_in_z(domain, dtypes + 1, array)){
      return 1;
    }
  }
  return 0;
}

int compute_flux_z(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict uz = fluid->uz.data;
  const double * restrict vof = interface->vof.data;
  const normal_t * restrict normal = interface->normal.data;
  double * restrict fluxz = interface->fluxz.data;
  for(int k = 1; k <= ksize + 1; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // use upwind information
        const double vel = UZ(i, j, k);
        const int    kk = vel < 0. ?    k : k - 1;
        const double  z = vel < 0. ? -0.5 :  +0.5;
        // evaluate flux
        const double lvof = VOF(i, j, kk);
        if(lvof < vofmin || 1. - vofmin < lvof){
          FLUXZ(i, j, k) = vel * lvof;
          continue;
        }
        double flux = 0.;
        for(int jj = 0; jj < NGAUSS; jj++){
          for(int ii = 0; ii < NGAUSS; ii++){
            const double w = gauss_ws[ii] * gauss_ws[jj];
            const double x = gauss_ps[ii];
            const double y = gauss_ps[jj];
            flux += w * indicator(NORMAL(i, j, kk), (const double [NDIMS]){x, y, z});
          }
        }
        FLUXZ(i, j, k) = vel * flux;
      }
    }
  }
  update_boundaries(domain, &interface->fluxz);
  return 0;
}
