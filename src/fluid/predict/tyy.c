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
  const double hy = domain->hy;
  const double * restrict uy = fluid->uy.data;
  const double * restrict visuy = fluid->visuy.data;
  double * restrict tyy = fluid->tyy.data;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      // tyy at cell center
      const double vis_uy = + 0.5 * VISUY(i  , j  )
                            + 0.5 * VISUY(i  , j+1);
      const double duy = - UY(i  , j  )
                         + UY(i  , j+1);
      const double lyy = 1. / hy * duy;
      TYY(i, j) = vis_uy * lyy + vis_uy * lyy;
    }
  }
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, &fluid->tyy)){
    return 1;
  }
  return 0;
}

