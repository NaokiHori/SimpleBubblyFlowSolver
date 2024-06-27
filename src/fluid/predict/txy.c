#include <mpi.h>
#include "domain.h"
#include "fluid.h"
#include "halo.h"
#include "./internal.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/visux.h"
#include "array_macros/fluid/visuy.h"
#include "array_macros/fluid/txy.h"

int compute_txy(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double hy = domain->hy;
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
  const double * restrict visux = fluid->visux.data;
  const double * restrict visuy = fluid->visuy.data;
  double * restrict txy = fluid->txy.data;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize + 1; i++){
      // txy at cell corner
      const double vis_ux = + 0.5 * VISUX(i  , j-1)
                            + 0.5 * VISUX(i  , j  );
      const double vis_uy = + 0.5 * VISUY(i-1, j  )
                            + 0.5 * VISUY(i  , j  );
      const double hx = HXXF(i  );
      const double dux = - UX(i  , j-1)
                         + UX(i  , j  );
      const double duy = - UY(i-1, j  )
                         + UY(i  , j  );
      const double lxy = 1. / hy * dux;
      const double lyx = 1. / hx * duy;
      TXY(i, j) = vis_ux * lxy + vis_uy * lyx;
    }
  }
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, &fluid->txy)){
    return 1;
  }
  return 0;
}

