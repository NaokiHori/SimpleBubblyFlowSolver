#include <math.h>
#include <float.h>
#include "param.h"
#include "domain.h"
#include "interface.h"
#include "internal.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/hxxc.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/dvof.h"
#include "array_macros/interface/normal.h"
#include "array_macros/interface/curv.h"

static int compute_gradient(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
  const double hy = domain->hy;
#if NDIMS == 3
  const double hz = domain->hz;
#endif
  const double * restrict vof = interface->vof.data;
  vector_t * restrict dvof = interface->dvof.data;
#if NDIMS == 2
  for(int j = 0; j <= jsize + 2; j++){
    for(int i = 1; i <= isize + 1; i++){
      // x gradient
      const double gradx = 1. / HXXF(i  ) * (
          - VOF(i-1, j-1) + VOF(i  , j-1)
          - VOF(i-1, j  ) + VOF(i  , j  )
      );
      // y gradient
      const double grady = 1. / hy * (
          - VOF(i-1, j-1) - VOF(i  , j-1)
          + VOF(i-1, j  ) + VOF(i  , j  )
      );
      // normalise and obtain corner normals
      const double norm = sqrt(
          + pow(gradx, 2.)
          + pow(grady, 2.)
      );
      const double norminv = 1. / fmax(norm, DBL_EPSILON);
      DVOF(i, j)[0] = gradx * norminv;
      DVOF(i, j)[1] = grady * norminv;
    }
  }
#else
  for(int k = 0; k <= ksize + 2; k++){
    for(int j = 0; j <= jsize + 2; j++){
      for(int i = 1; i <= isize + 1; i++){
        // x gradient
        const double gradx = 1. / HXXF(i  ) * (
            - VOF(i-1, j-1, k-1) + VOF(i  , j-1, k-1)
            - VOF(i-1, j  , k-1) + VOF(i  , j  , k-1)
            - VOF(i-1, j-1, k  ) + VOF(i  , j-1, k  )
            - VOF(i-1, j  , k  ) + VOF(i  , j  , k  )
        );
        // y gradient
        const double grady = 1. / hy * (
            - VOF(i-1, j-1, k-1) - VOF(i  , j-1, k-1)
            + VOF(i-1, j  , k-1) + VOF(i  , j  , k-1)
            - VOF(i-1, j-1, k  ) - VOF(i  , j-1, k  )
            + VOF(i-1, j  , k  ) + VOF(i  , j  , k  )
        );
        // z gradient
        const double gradz = 1. / hz * (
            - VOF(i-1, j-1, k-1) - VOF(i  , j-1, k-1)
            - VOF(i-1, j  , k-1) - VOF(i  , j  , k-1)
            + VOF(i-1, j-1, k  ) + VOF(i  , j-1, k  )
            + VOF(i-1, j  , k  ) + VOF(i  , j  , k  )
        );
        // normalise and obtain corner normals
        const double norm = sqrt(
            + pow(gradx, 2.)
            + pow(grady, 2.)
            + pow(gradz, 2.)
        );
        const double norminv = 1. / fmax(norm, DBL_EPSILON);
        DVOF(i, j, k)[0] = gradx * norminv;
        DVOF(i, j, k)[1] = grady * norminv;
        DVOF(i, j, k)[2] = gradz * norminv;
      }
    }
  }
#endif
  return 0;
}

static int compute_normal(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxc = domain->hxxc;
  const double hy = domain->hy;
#if NDIMS == 3
  const double hz = domain->hz;
#endif
  const double * restrict vof = interface->vof.data;
  const vector_t * restrict dvof = interface->dvof.data;
  normal_t * restrict normal = interface->normal.data;
#if NDIMS == 2
  for(int j = 0; j <= jsize + 1; j++){
    for(int i = 1; i <= isize; i++){
      const double lvof = VOF(i, j);
      // for (almost) single-phase region,
      //   surface reconstruction is not needed
      if(lvof < vofmin || 1. - vofmin < lvof){
        continue;
      }
      // average nx
      double nx = (
          + DVOF(i  , j  )[0] + DVOF(i+1, j  )[0]
          + DVOF(i  , j+1)[0] + DVOF(i+1, j+1)[0]
      );
      // average ny
      double ny = (
          + DVOF(i  , j  )[1] + DVOF(i+1, j  )[1]
          + DVOF(i  , j+1)[1] + DVOF(i+1, j+1)[1]
      );
      // normalise and obtain center normals
      nx /= HXXC(i  );
      ny /= hy;
      const double norm = sqrt(
          + pow(nx, 2.)
          + pow(ny, 2.)
      );
      const double norminv = 1. / fmax(norm, DBL_EPSILON);
      nx *= norminv;
      ny *= norminv;
      // store normal and intercept
      NORMAL(i, j)[0] = nx;
      NORMAL(i, j)[1] = ny;
      NORMAL(i, j)[2] = - 0.5 / vofbeta * log(1. / lvof - 1.);
    }
  }
#else
  for(int k = 0; k <= ksize + 1; k++){
    for(int j = 0; j <= jsize + 1; j++){
      for(int i = 1; i <= isize; i++){
        const double lvof = VOF(i, j, k);
        // for (almost) single-phase region,
        //   surface reconstruction is not needed
        if(lvof < vofmin || 1. - vofmin < lvof){
          continue;
        }
        // average nx
        double nx = (
            + DVOF(i  , j  , k  )[0] + DVOF(i+1, j  , k  )[0]
            + DVOF(i  , j+1, k  )[0] + DVOF(i+1, j+1, k  )[0]
            + DVOF(i  , j  , k+1)[0] + DVOF(i+1, j  , k+1)[0]
            + DVOF(i  , j+1, k+1)[0] + DVOF(i+1, j+1, k+1)[0]
        );
        // average ny
        double ny = (
            + DVOF(i  , j  , k  )[1] + DVOF(i+1, j  , k  )[1]
            + DVOF(i  , j+1, k  )[1] + DVOF(i+1, j+1, k  )[1]
            + DVOF(i  , j  , k+1)[1] + DVOF(i+1, j  , k+1)[1]
            + DVOF(i  , j+1, k+1)[1] + DVOF(i+1, j+1, k+1)[1]
        );
        // average nz
        double nz = (
            + DVOF(i  , j  , k  )[2] + DVOF(i+1, j  , k  )[2]
            + DVOF(i  , j+1, k  )[2] + DVOF(i+1, j+1, k  )[2]
            + DVOF(i  , j  , k+1)[2] + DVOF(i+1, j  , k+1)[2]
            + DVOF(i  , j+1, k+1)[2] + DVOF(i+1, j+1, k+1)[2]
        );
        // normalise and obtain center normals
        nx = nx / HXXC(i  );
        ny = ny / hy;
        nz = nz / hz;
        const double norm = sqrt(
            + pow(nx, 2.)
            + pow(ny, 2.)
            + pow(nz, 2.)
        );
        const double norminv = 1. / fmax(norm, DBL_EPSILON);
        nx *= norminv;
        ny *= norminv;
        nz *= norminv;
        // store normal and intercept
        NORMAL(i, j, k)[0] = nx;
        NORMAL(i, j, k)[1] = ny;
        NORMAL(i, j, k)[2] = nz;
        NORMAL(i, j, k)[3] = - 0.5 / vofbeta * log(1. / lvof - 1.);
      }
    }
  }
#endif
  return 0;
}

// compute mean curvature from corner normals
static int compute_curvature(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
  const double hy = domain->hy;
#if NDIMS == 3
  const double hz = domain->hz;
#endif
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  const vector_t * restrict dvof = interface->dvof.data;
  double * restrict curv = interface->curv.data;
#if NDIMS == 2
  for(int j = 0; j <= jsize + 1; j++){
    for(int i = 1; i <= isize; i++){
      const double hx_xm = HXXF(i  );
      const double hx_xp = HXXF(i+1);
      const double jd_xm = JDXF(i  );
      const double jd_x0 = JDXC(i  );
      const double jd_xp = JDXF(i+1);
      const double nx_xm = + 0.5 * DVOF(i  , j  )[0]
                           + 0.5 * DVOF(i  , j+1)[0];
      const double nx_xp = + 0.5 * DVOF(i+1, j  )[0]
                           + 0.5 * DVOF(i+1, j+1)[0];
      const double ny_ym = + 0.5 * DVOF(i  , j  )[1]
                           + 0.5 * DVOF(i+1, j  )[1];
      const double ny_yp = + 0.5 * DVOF(i  , j+1)[1]
                           + 0.5 * DVOF(i+1, j+1)[1];
      const double div = 1. / jd_x0 * (
          - jd_xm / hx_xm * nx_xm + jd_xp / hx_xp * nx_xp
          - jd_x0 / hy    * ny_ym + jd_x0 / hy    * ny_yp
      );
      CURV(i, j) = - 1. * div;
    }
  }
#else
  for(int k = 0; k <= ksize + 1; k++){
    for(int j = 0; j <= jsize + 1; j++){
      for(int i = 1; i <= isize; i++){
        const double hx_xm = HXXF(i  );
        const double hx_xp = HXXF(i+1);
        const double jd_xm = JDXF(i  );
        const double jd_x0 = JDXC(i  );
        const double jd_xp = JDXF(i+1);
        const double nx_xm = + 0.25 * DVOF(i  , j  , k  )[0]
                             + 0.25 * DVOF(i  , j+1, k  )[0]
                             + 0.25 * DVOF(i  , j  , k+1)[0]
                             + 0.25 * DVOF(i  , j+1, k+1)[0];
        const double nx_xp = + 0.25 * DVOF(i+1, j  , k  )[0]
                             + 0.25 * DVOF(i+1, j+1, k  )[0]
                             + 0.25 * DVOF(i+1, j  , k+1)[0]
                             + 0.25 * DVOF(i+1, j+1, k+1)[0];
        const double ny_ym = + 0.25 * DVOF(i  , j  , k  )[1]
                             + 0.25 * DVOF(i+1, j  , k  )[1]
                             + 0.25 * DVOF(i  , j  , k+1)[1]
                             + 0.25 * DVOF(i+1, j  , k+1)[1];
        const double ny_yp = + 0.25 * DVOF(i  , j+1, k  )[1]
                             + 0.25 * DVOF(i+1, j+1, k  )[1]
                             + 0.25 * DVOF(i  , j+1, k+1)[1]
                             + 0.25 * DVOF(i+1, j+1, k+1)[1];
        const double nz_zm = + 0.25 * DVOF(i  , j  , k  )[2]
                             + 0.25 * DVOF(i+1, j  , k  )[2]
                             + 0.25 * DVOF(i  , j+1, k  )[2]
                             + 0.25 * DVOF(i+1, j+1, k  )[2];
        const double nz_zp = + 0.25 * DVOF(i  , j  , k+1)[2]
                             + 0.25 * DVOF(i+1, j  , k+1)[2]
                             + 0.25 * DVOF(i  , j+1, k+1)[2]
                             + 0.25 * DVOF(i+1, j+1, k+1)[2];
        const double div = 1. / jd_x0 * (
            - jd_xm / hx_xm * nx_xm + jd_xp / hx_xp * nx_xp
            - jd_x0 / hy    * ny_ym + jd_x0 / hy    * ny_yp
            - jd_x0 / hz    * nz_zm + jd_x0 / hz    * nz_zp
        );
        CURV(i, j, k) = - 1. * div;
      }
    }
  }
#endif
  return 0;
}

int interface_compute_curvature_tensor(
    const domain_t * domain,
    interface_t * interface
){
  compute_gradient(domain, interface);
  compute_normal(domain, interface);
  compute_curvature(domain, interface);
  return 0;
}

