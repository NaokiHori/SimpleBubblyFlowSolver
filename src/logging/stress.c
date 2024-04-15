#include <stdio.h>
#include <math.h>
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/fluid/txy.h"
#include "internal.h"

int logging_check_stress(
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
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict txy = fluid->txy.data;
  const double diffusivity = 1. / fluid->Re;
  // shear stress in the y direction for each wall
  double vals[2] = {0., 0.};
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    const int im =         1;
    const int ip = isize + 1;
    vals[0] += diffusivity * JDXF(im) / HXXF(im) * TXY(im, j);
    vals[1] += diffusivity * JDXF(ip) / HXXF(ip) * TXY(ip, j);
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      const int im =         1;
      const int ip = isize + 1;
      vals[0] += diffusivity * JDXF(im) / HXXF(im) * TXY(im, j, k);
      vals[1] += diffusivity * JDXF(ip) / HXXF(ip) * TXY(ip, j, k);
    }
  }
#endif
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : vals;
  void * recvbuf = vals;
  MPI_Reduce(sendbuf, recvbuf, 2, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  if(root == myrank){
    FILE * fp = fileio.fopen(fname, "a");
    if(NULL == fp){
      return 0;
    }
    fprintf(fp, "%8.2f % 18.15e % 18.15e\n", time, vals[0], vals[1]);
    fileio.fclose(fp);
  }
  return 0;
}

