#include <stdio.h>
#include <math.h>
#include "domain.h"
#include "interface.h"
#include "fileio.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/interface/vof.h"
#include "internal.h"

int logging_check_vof(
    const char fname[],
    const domain_t * domain,
    const double time,
    const interface_t * interface
){
  const int root = 0;
  int myrank = root;
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_rank(domain->info, &myrank);
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * jdxf = domain->jdxf;
  const double * jdxc = domain->jdxc;
  const double * hxxf = domain->hxxf;
  const double hy = domain->hy;
  const double hz = domain->hz;
  const double * vof = interface->vof.data;
  double min = 1.;
  double max = 0.;
  double sums[3] = {0.};
  // check volume conservation
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double dv = JDXC(i  );
        const double lvof = VOF(i, j, k);
        min = fmin(min, lvof);
        max = fmax(max, lvof);
        sums[0] += lvof * dv;
        sums[1] +=        dv;
      }
    }
  }
  // surface area, x face
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize + 1; i++){
        const double vof_xm = VOF(i-1, j  , k  );
        const double vof_xp = VOF(i  , j  , k  );
        const double ds = JDXF(i  ) / HXXF(i  );
        sums[2] += fabs(vof_xp - vof_xm) * ds;
      }
    }
  }
  // surface area, y face
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double vof_ym = VOF(i  , j-1, k  );
        const double vof_yp = VOF(i  , j  , k  );
        const double ds = JDXC(i  ) / hy;
        sums[2] += fabs(vof_yp - vof_ym) * ds;
      }
    }
  }
  // surface area, z face
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double vof_zm = VOF(i  , j  , k-1);
        const double vof_zp = VOF(i  , j  , k  );
        const double ds = JDXC(i  ) / hz;
        sums[2] += fabs(vof_zp - vof_zm) * ds;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DOUBLE, MPI_MIN, comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE, MPI_MAX, comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, sums, 3, MPI_DOUBLE, MPI_SUM, comm_cart);
  if(root == myrank){
    FILE * fp = fileio.fopen(fname, "a");
    if(NULL == fp){
      return 0;
    }
    fprintf(fp, "%8.2f ", time);
    fprintf(fp, "% 18.15e ", min);
    fprintf(fp, "% 18.15e ", max);
    fprintf(fp, "% 18.15e ", sums[0] / sums[1]);
    fprintf(fp, "% 18.15e\n", sums[2]);
    fileio.fclose(fp);
  }
  return 0;
}

