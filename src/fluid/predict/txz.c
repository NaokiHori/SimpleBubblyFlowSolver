#if NDIMS == 3
#include <mpi.h>
#include "domain.h"
#include "fluid.h"
#include "halo.h"
#include "./internal.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/visux.h"
#include "array_macros/fluid/visuz.h"
#include "array_macros/fluid/txz.h"

int compute_txz(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxf = domain->hxxf;
  const double hz = domain->hz;
  const double * restrict ux = fluid->ux.data;
  const double * restrict uz = fluid->uz.data;
  const double * restrict visux = fluid->visux.data;
  const double * restrict visuz = fluid->visuz.data;
  double * restrict txz = fluid->txz.data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize + 1; i++){
        // txz at cell corner | 12
        const double vis_ux = + 0.5 * VISUX(i  , j  , k-1)
                              + 0.5 * VISUX(i  , j  , k  );
        const double vis_uz = + 0.5 * VISUZ(i-1, j  , k  )
                              + 0.5 * VISUZ(i  , j  , k  );
        const double hx = HXXF(i  );
        const double dux = - UX(i  , j  , k-1)
                           + UX(i  , j  , k  );
        const double duz = - UZ(i-1, j  , k  )
                           + UZ(i  , j  , k  );
        const double lxz = 1. / hz * dux;
        const double lzx = 1. / hx * duz;
        TXZ(i, j, k) = vis_ux * lxz + vis_uz * lzx;
      }
    }
  }
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, &fluid->txz)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, &fluid->txz)){
    return 1;
  }
  return 0;
}
#endif
