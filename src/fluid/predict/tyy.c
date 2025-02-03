#include <mpi.h>
#include "domain.h"
#include "fluid.h"
#include "halo.h"
#include "./internal.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/visuy.h"
#include "array_macros/fluid/tyy.h"

int compute_tyy(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hy = domain->hy;
  const double * restrict uy = fluid->uy.data;
  const double * restrict visuy = fluid->visuy.data;
  double * restrict tyy = fluid->tyy.data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // tyy at cell center
        const double vis_uy = + 0.5 * VISUY(i  , j  , k  )
                              + 0.5 * VISUY(i  , j+1, k  );
        const double duy = - UY(i  , j  , k  )
                           + UY(i  , j+1, k  );
        const double lyy = 1. / hy * duy;
        TYY(i, j, k) = vis_uy * lyy + vis_uy * lyy;
      }
    }
  }
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, &fluid->tyy)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, &fluid->tyy)){
    return 1;
  }
  return 0;
}

