#include <assert.h>
#include <mpi.h>
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "interface_solver.h"
#include "array_macros/interface/vof.h"

/**
 * @brief update boundary values of vof field
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] array  : volume-of-fluid
 * @return               : error code
 */
int interface_update_boundaries_vof(
    const domain_t * domain,
    array_t * array
){
  for (size_t dim = 0; dim < NDIMS; dim++) {
    assert(VOF_NADDS[dim][0] == array->nadds[dim][0]);
    assert(VOF_NADDS[dim][1] == array->nadds[dim][1]);
  }
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    double * vof = array->data;
    // set boundary values
    for(int j = 1; j <= jsize; j++){
      VOF(      0, j) = 0.;
      VOF(isize+1, j) = 0.;
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

