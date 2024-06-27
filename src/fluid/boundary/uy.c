#include <assert.h>
#include <mpi.h>
#include "param.h"
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/uy.h"

/**
 * @brief update boundary values of y velocity
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] array  : y velocity
 * @return               : error code
 */
int fluid_update_boundaries_uy(
    const domain_t * domain,
    array_t * array
){
  for (size_t dim = 0; dim < NDIMS; dim++) {
    assert(UY_NADDS[dim][0] == array->nadds[dim][0]);
    assert(UY_NADDS[dim][1] == array->nadds[dim][1]);
  }
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const int ksize = domain->mysizes[2];
    double * uy = array->data;
    // set boundary values
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        UY(      0, j, k) = param_uy_xm;
        UY(isize+1, j, k) = param_uy_xp;
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

