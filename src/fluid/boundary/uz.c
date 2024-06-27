#include <assert.h>
#include <mpi.h>
#include "param.h"
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/uz.h"

/**
 * @brief update boundary values of z velocity
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] array  : z velocity
 * @return               : error code
 */
int fluid_update_boundaries_uz(
    const domain_t * domain,
    array_t * array
){
  for (size_t dim = 0; dim < NDIMS; dim++) {
    assert(UZ_NADDS[dim][0] == array->nadds[dim][0]);
    assert(UZ_NADDS[dim][1] == array->nadds[dim][1]);
  }
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const int ksize = domain->mysizes[2];
    double * uz = array->data;
    // set boundary values
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        UZ(      0, j, k) = param_uz_xm;
        UZ(isize+1, j, k) = param_uz_xp;
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
