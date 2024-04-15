#if NDIMS == 3
#include <mpi.h>
#include "domain.h"
#include "fluid.h"
#include "halo.h"
#include "./internal.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/visuz.h"
#include "array_macros/fluid/tzz.h"

int compute_tzz(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict uz = fluid->uz.data;
  const double * restrict visuz = fluid->visuz.data;
  double * restrict tzz = fluid->tzz.data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // tzz at cell center | 6
        const double vis_uz = + 0.5 * VISUZ(i  , j  , k  )
                              + 0.5 * VISUZ(i  , j  , k+1);
        const double duz = - UZ(i  , j  , k  )
                           + UZ(i  , j  , k+1);
        const double lzz = 1. / hz * duz;
        TZZ(i, j, k) = vis_uz * lzz + vis_uz * lzz;
      }
    }
  }
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, &fluid->tzz)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, &fluid->tzz)){
    return 1;
  }
  return 0;
}
#endif
