#include <mpi.h>
#include "domain.h"
#include "fluid.h"
#include "halo.h"
#include "./internal.h"
#include "array_macros/domain/hxxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/visux.h"
#include "array_macros/fluid/txx.h"

int compute_txx(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxc = domain->hxxc;
  const double * restrict ux = fluid->ux.data;
  const double * restrict visux = fluid->visux.data;
  double * restrict txx = fluid->txx.data;
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      // txx at cell center | 7
      const double vis_ux = + 0.5 * VISUX(i  , j  )
                            + 0.5 * VISUX(i+1, j  );
      const double hx = HXXC(i  );
      const double dux = - UX(i  , j  )
                         + UX(i+1, j  );
      const double lxx = 1. / hx * dux;
      TXX(i, j) = vis_ux * lxx + vis_ux * lxx;
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // txx at cell center | 7
        const double vis_ux = + 0.5 * VISUX(i  , j  , k  )
                              + 0.5 * VISUX(i+1, j  , k  );
        const double hx = HXXC(i  );
        const double dux = - UX(i  , j  , k  )
                           + UX(i+1, j  , k  );
        const double lxx = 1. / hx * dux;
        TXX(i, j, k) = vis_ux * lxx + vis_ux * lxx;
      }
    }
  }
#endif
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
#if NDIMS == 3
    MPI_DOUBLE,
#endif
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, &fluid->txx)){
    return 1;
  }
#if NDIMS == 3
  if(0 != halo_communicate_in_z(domain, dtypes + 1, &fluid->txx)){
    return 1;
  }
#endif
  return 0;
}

