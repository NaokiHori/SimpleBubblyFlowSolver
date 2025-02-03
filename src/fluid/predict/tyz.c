#include <mpi.h>
#include "domain.h"
#include "fluid.h"
#include "halo.h"
#include "./internal.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/visuy.h"
#include "array_macros/fluid/visuz.h"
#include "array_macros/fluid/tyz.h"

int compute_tyz(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hy = domain->hy;
  const double hz = domain->hz;
  const double * restrict uy = fluid->uy.data;
  const double * restrict uz = fluid->uz.data;
  const double * restrict visuy = fluid->visuy.data;
  const double * restrict visuz = fluid->visuz.data;
  double * restrict tyz = fluid->tyz.data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // tyz at cell corner
        const double vis_uy = + 0.5 * VISUY(i  , j  , k-1)
                              + 0.5 * VISUY(i  , j  , k  );
        const double vis_uz = + 0.5 * VISUZ(i  , j-1, k  )
                              + 0.5 * VISUZ(i  , j  , k  );
        const double duy = - UY(i  , j  , k-1)
                           + UY(i  , j  , k  );
        const double duz = - UZ(i  , j-1, k  )
                           + UZ(i  , j  , k  );
        const double lyz = 1. / hz * duy;
        const double lzy = 1. / hy * duz;
        TYZ(i, j, k) = vis_uy * lyz + vis_uz * lzy;
      }
    }
  }
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, &fluid->tyz)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, &fluid->tyz)){
    return 1;
  }
  return 0;
}
