#if NDIMS == 3
#include "memory.h"
#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "interface.h"
#include "interface_solver.h"
#include "internal.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/p.h"
#include "array_macros/fluid/den.h"
#include "array_macros/fluid/txz.h"
#include "array_macros/fluid/tyz.h"
#include "array_macros/fluid/tzz.h"
#include "array_macros/interface/fluxx.h"
#include "array_macros/interface/fluxy.h"
#include "array_macros/interface/fluxz.h"
#include "array_macros/interface/ifrcz.h"

#define BEGIN \
  for(int cnt = 0, k = 1; k <= ksize; k++){ \
    for(int j = 1; j <= jsize; j++){ \
      for(int i = 1; i <= isize; i++, cnt++){
#define END \
      } \
    } \
  }

static int advection_x(
    const domain_t * domain,
    const double * restrict uz,
    const double * restrict fluxx,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    // uz is advected in x | 17
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
    const double mux_xm = + 0.5 * jd_xm / hx_xm * FLUXX(i  , j  , k-1)
                          + 0.5 * jd_xm / hx_xm * FLUXX(i  , j  , k  );
    const double mux_xp = + 0.5 * jd_xp / hx_xp * FLUXX(i+1, j  , k-1)
                          + 0.5 * jd_xp / hx_xp * FLUXX(i+1, j  , k  );
    const double uz_xm = + 0.5 * UZ(i-1, j  , k  )
                         + 0.5 * UZ(i  , j  , k  );
    const double uz_xp = + 0.5 * UZ(i  , j  , k  )
                         + 0.5 * UZ(i+1, j  , k  );
    src[cnt] -= 1. / jd_x0 * (
        - mux_xm * uz_xm
        + mux_xp * uz_xp
    );
  END
  return 0;
}

static int advection_y(
    const domain_t * domain,
    const double * restrict uz,
    const double * restrict fluxy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hy = domain->hy;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    // uz is advected in y | 13
    const double jd = JDXC(i  );
    const double muy_ym = + 0.5 * jd / hy * FLUXY(i  , j  , k-1)
                          + 0.5 * jd / hy * FLUXY(i  , j  , k  );
    const double muy_yp = + 0.5 * jd / hy * FLUXY(i  , j+1, k-1)
                          + 0.5 * jd / hy * FLUXY(i  , j+1, k  );
    const double uz_ym = + 0.5 * UZ(i  , j-1, k  )
                         + 0.5 * UZ(i  , j  , k  );
    const double uz_yp = + 0.5 * UZ(i  , j  , k  )
                         + 0.5 * UZ(i  , j+1, k  );
    src[cnt] -= 1. / jd * (
        - muy_ym * uz_ym
        + muy_yp * uz_yp
    );
  END
  return 0;
}

static int advection_z(
    const domain_t * domain,
    const double * restrict uz,
    const double * restrict fluxz,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    // uz is advected in z | 13
    const double jd = JDXC(i  );
    const double muz_zm = + 0.5 * jd / hz * FLUXZ(i  , j  , k-1)
                          + 0.5 * jd / hz * FLUXZ(i  , j  , k  );
    const double muz_zp = + 0.5 * jd / hz * FLUXZ(i  , j  , k  )
                          + 0.5 * jd / hz * FLUXZ(i  , j  , k+1);
    const double uz_zm = + 0.5 * UZ(i  , j  , k-1)
                         + 0.5 * UZ(i  , j  , k  );
    const double uz_zp = + 0.5 * UZ(i  , j  , k  )
                         + 0.5 * UZ(i  , j  , k+1);
    src[cnt] -= 1. / jd * (
        - muz_zm * uz_zm
        + muz_zp * uz_zp
    );
  END
  return 0;
}

static int diffusion_x(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict txz,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    // uz is diffused in x | 11
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
    const double txz_xm = TXZ(i  , j  , k  );
    const double txz_xp = TXZ(i+1, j  , k  );
    src[cnt] += diffusivity / jd_x0 * (
        - jd_xm / hx_xm * txz_xm
        + jd_xp / hx_xp * txz_xp
    );
  END
  return 0;
}

static int diffusion_y(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict tyz,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hy = domain->hy;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    // uz is diffused in y | 7
    const double jd = JDXC(i  );
    const double tyz_ym = TYZ(i  , j  , k  );
    const double tyz_yp = TYZ(i  , j+1, k  );
    src[cnt] += diffusivity / jd * (
        - jd / hy * tyz_ym
        + jd / hy * tyz_yp
    );
  END
  return 0;
}

static int diffusion_z(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict tzz,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    // uz is diffused in z | 7
    const double jd = JDXC(i  );
    const double tzz_zm = TZZ(i  , j  , k-1);
    const double tzz_zp = TZZ(i  , j  , k  );
    src[cnt] += diffusivity / jd * (
        - jd / hz * tzz_zm
        + jd / hz * tzz_zp
    );
  END
  return 0;
}

static int pressure(
    const domain_t * domain,
    const double * restrict p,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  BEGIN
    // pressure-gradient contribution | 4
    src[cnt] -= 1. / hz * (
        - P(i  , j  , k-1)
        + P(i  , j  , k  )
    );
  END
  return 0;
}

static int surface(
    const domain_t * domain,
    const double * restrict ifrcz,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  BEGIN
    src[cnt] += IFRCZ(i, j, k);
  END
  return 0;
}

/**
 * @brief comute right-hand-side of Runge-Kutta scheme of uz
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in,out] fluid  : n-step flow field (in), RK source terms (inout)
 * @return               : error code
 */
int compute_rhs_uz(
    const domain_t * domain,
    fluid_t * fluid,
    const interface_t * interface
){
  const double * restrict  uz = fluid-> uz.data;
  const double * restrict   p = fluid->  p.data;
  const double * restrict txz = fluid->txz.data;
  const double * restrict tyz = fluid->tyz.data;
  const double * restrict tzz = fluid->tzz.data;
  const double * restrict fluxx = interface->fluxx.data;
  const double * restrict fluxy = interface->fluxy.data;
  const double * restrict fluxz = interface->fluxz.data;
  double * restrict srca = fluid->srcuz[rk_a].data;
  double * restrict srcg = fluid->srcuz[rk_g].data;
  const double diffusivity = 1. / fluid->Re;
  // advective contributions
  advection_x(domain, uz, fluxx, srca);
  advection_y(domain, uz, fluxy, srca);
  advection_z(domain, uz, fluxz, srca);
  // diffusive contributions
  diffusion_x(domain, diffusivity, txz, srca);
  diffusion_y(domain, diffusivity, tyz, srca);
  diffusion_z(domain, diffusivity, tzz, srca);
  // pressure-gradient contribution
  pressure(domain, p, srcg);
  surface(domain, interface->ifrcz.data, srca);
  return 0;
}

/**
 * @brief predict uz
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : Runge-Kutta source terms (in), velocity (out)
 * @return               : error code
 */
int predict_uz(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double coef_a = rkcoefs[rkstep][rk_a];
  const double coef_b = rkcoefs[rkstep][rk_b];
  const double coef_g = rkcoefs[rkstep][rk_g];
  const double * restrict srcuza = fluid->srcuz[rk_a].data;
  const double * restrict srcuzb = fluid->srcuz[rk_b].data;
  const double * restrict srcuzg = fluid->srcuz[rk_g].data;
  double * restrict uz = fluid->uz.data;
  {
    const double * restrict den = fluid->den[0].data;
    BEGIN
      const double lden = + 0.5 * DEN(i  , j  , k-1)
                          + 0.5 * DEN(i  , j  , k  );
      UZ(i, j, k) =
        + lden * UZ(i, j, k)
        + coef_a * dt * srcuza[cnt]
        + coef_b * dt * srcuzb[cnt]
        + coef_g * dt * srcuzg[cnt];
    END
  }
  {
    const double * restrict den = fluid->den[1].data;
    BEGIN
      const double lden = + 0.5 * DEN(i  , j  , k-1)
                          + 0.5 * DEN(i  , j  , k  );
      UZ(i, j, k) /= lden;
    END
  }
  fluid_update_boundaries_uz(domain, &fluid->uz);
  return 0;
}
#endif
