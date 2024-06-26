#include <math.h>
#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#include "interface.h"
#include "interface_solver.h"
#include "../internal.h"
#include "internal.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/fluxx.h"
#include "array_macros/interface/fluxy.h"
#include "array_macros/interface/src.h"

// planar surface
inline static double surface_function(
    const normal_t n,
    const vector_t x
){
  return
    + n[0] * x[0]
    + n[1] * x[1]
    + n[NDIMS];
}

// diffused surface representation
double indicator(
    const normal_t n,
    const vector_t x
){
  const double sf = surface_function(n, x);
  return 0.5 * (1. + tanh(vofbeta * sf));
}

static int reset_srcs(
    const size_t rkstep,
    array_t * restrict srca,
    array_t * restrict srcb
){
  // copy previous k-step source term and reset
  if(0 != rkstep){
    // stash previous RK source term,
    //   which is achieved by swapping
    //   the pointers to "data"
    double * tmp = srca->data;
    srca->data = srcb->data;
    srcb->data = tmp;
  }
  return 0;
}

static int interface_compute_rhs(
    const domain_t * domain,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double hy = domain->hy;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict fluxx = interface->fluxx.data;
  const double * restrict fluxy = interface->fluxy.data;
  double * restrict src = interface->src[rk_a].data;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      // compute source of volume-of-fluid
      const double hx_xm = HXXF(i  );
      const double hx_xp = HXXF(i+1);
      const double jd_xm = JDXF(i  );
      const double jd_x0 = JDXC(i  );
      const double jd_xp = JDXF(i+1);
      const double flux_xm = FLUXX(i  , j  );
      const double flux_xp = FLUXX(i+1, j  );
      const double flux_ym = FLUXY(i  , j  );
      const double flux_yp = FLUXY(i  , j+1);
      SRC(i, j) = 1. / jd_x0 * (
          + jd_xm / hx_xm * flux_xm
          - jd_xp / hx_xp * flux_xp
          + jd_x0 / hy    * flux_ym
          - jd_x0 / hy    * flux_yp
      );
    }
  }
  return 0;
}

static int interface_advect_vof(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    interface_t * interface
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  double * restrict vof = interface->vof.data;
  // update vof, alpha contribution
  {
    const double coef = rkcoefs[rkstep][rk_a];
    const double * restrict src = interface->src[rk_a].data;
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        VOF(i, j) += dt * coef * SRC(i, j);
      }
    }
  }
  // update vof, beta contribution
  if(0 != rkstep){
    const double coef = rkcoefs[rkstep][rk_b];
    const double * restrict src = interface->src[rk_b].data;
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        VOF(i, j) += dt * coef * SRC(i, j);
      }
    }
  }
  return 0;
}

int interface_update_vof(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    const fluid_t * fluid,
    interface_t * interface
){
  reset_srcs(rkstep, interface->src + rk_a, interface->src + rk_b);
  compute_flux_x(domain, fluid, interface);
  compute_flux_y(domain, fluid, interface);
  interface_compute_rhs(domain, interface);
  interface_advect_vof(domain, rkstep, dt, interface);
  interface_update_boundaries_vof(domain, &interface->vof);
  return 0;
}

