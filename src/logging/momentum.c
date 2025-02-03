#include <stdio.h>
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/den.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "internal.h"

/**
 * @brief compute total momenta
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information about domain decomposition and size
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int logging_check_momentum(
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
){
  const int root = 0;
  int myrank = root;
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_rank(domain->info, &myrank);
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict den = fluid->den[1].data;
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
  double moms[NDIMS] = {0.};
  // compute total x-momentum
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      const double ds = JDXF(i  );
      const double lden = 0.5 * DEN(i-1, j  ) + 0.5 * DEN(i  , j  );
      const double lvel = UX(i, j);
      moms[0] += lden * lvel * ds;
    }
  }
  // compute total y-momentum
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      const double ds = JDXC(i  );
      const double lden = 0.5 * DEN(i  , j-1) + 0.5 * DEN(i  , j  );
      const double lvel = UY(i, j);
      moms[1] += lden * lvel * ds;
    }
  }
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : moms;
  void * recvbuf = moms;
  MPI_Reduce(sendbuf, recvbuf, NDIMS, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  if(root == myrank){
    FILE * fp = fileio.fopen(fname, "a");
    if(NULL == fp){
      return 0;
    }
    fprintf(fp, "%8.2f ", time);
    for(int n = 0; n < NDIMS; n++){
      fprintf(fp, "% 18.15e%c", moms[n], NDIMS - 1 == n ? '\n' : ' ');
    }
    fileio.fclose(fp);
  }
  return 0;
}

