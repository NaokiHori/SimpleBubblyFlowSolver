#include <assert.h>
#include <mpi.h>
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/p.h"

/**
 * @brief update boundary values of the pressure
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] array  : pressure
 * @return               : error code
 */
int fluid_update_boundaries_p(
    const domain_t * domain,
    array_t * array
){
  for (size_t dim = 0; dim < NDIMS; dim++) {
    assert(P_NADDS[dim][0] == array->nadds[dim][0]);
    assert(P_NADDS[dim][1] == array->nadds[dim][1]);
  }
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
#if NDIMS == 3
    const int ksize = domain->mysizes[2];
#endif
    double * p = array->data;
    // Neumann
#if NDIMS == 2
    for(int j = 1; j <= jsize; j++){
      P(      0, j) = P(    1, j);
      P(isize+1, j) = P(isize, j);
    }
#else
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        P(      0, j, k) = P(    1, j, k);
        P(isize+1, j, k) = P(isize, j, k);
      }
    }
#endif
  }
  {
    static MPI_Datatype dtypes[NDIMS - 1] = {
      MPI_DOUBLE,
#if NDIMS == 3
      MPI_DOUBLE,
#endif
    };
    if(0 != halo_communicate_in_y(domain, dtypes + 0, array)){
      return 1;
    }
#if NDIMS == 3
    if(0 != halo_communicate_in_z(domain, dtypes + 1, array)){
      return 1;
    }
#endif
  }
  return 0;
}

