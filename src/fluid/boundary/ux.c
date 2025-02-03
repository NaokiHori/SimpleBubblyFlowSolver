#include <assert.h>
#include <mpi.h>
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/ux.h"

/**
 * @brief update boundary values of x velocity
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] array  : x velocity
 * @return               : error code
 */
int fluid_update_boundaries_ux(
    const domain_t * domain,
    array_t * array
){
  for (size_t dim = 0; dim < NDIMS; dim++) {
    assert(UX_NADDS[dim][0] == array->nadds[dim][0]);
    assert(UX_NADDS[dim][1] == array->nadds[dim][1]);
  }
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    double * ux = array->data;
    // impermeable
    for(int j = 1; j <= jsize; j++){
      UX(      1, j) = 0.;
      UX(isize+1, j) = 0.;
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

