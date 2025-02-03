#include <stdio.h>
#include <math.h>
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/uz.h"
#include "internal.h"

int logging_check_extrema(
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
  const int ksize = domain->mysizes[2];
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
  const double * restrict uz = fluid->uz.data;
  double extrema[NDIMS] = {0.};
  // compute total x-momentum
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        extrema[0] = fmax(extrema[0], fabs(UX(i, j, k)));
      }
    }
  }
  // compute total y-momentum
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        extrema[1] = fmax(extrema[1], fabs(UY(i, j, k)));
      }
    }
  }
  // compute total z-momentum
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        extrema[2] = fmax(extrema[2], fabs(UZ(i, j, k)));
      }
    }
  }
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : extrema;
  void * recvbuf = extrema;
  MPI_Reduce(sendbuf, recvbuf, NDIMS, MPI_DOUBLE, MPI_MAX, root, comm_cart);
  if(root == myrank){
    FILE * fp = fileio.fopen(fname, "a");
    if(NULL == fp){
      return 0;
    }
    fprintf(fp, "%8.2f ", time);
    for(int n = 0; n < NDIMS; n++){
      fprintf(fp, "% 18.15e%c", extrema[n], NDIMS - 1 == n ? '\n' : ' ');
    }
    fileio.fclose(fp);
  }
  return 0;
}

