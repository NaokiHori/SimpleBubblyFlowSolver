#include "memory.h"
#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "interface.h"
#include "interface_solver.h"
#include "internal.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/hxxc.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/p.h"
#include "array_macros/fluid/den.h"
#include "array_macros/fluid/txx.h"
#include "array_macros/fluid/txy.h"
#if NDIMS == 3
#include "array_macros/fluid/txz.h"
#endif
#include "array_macros/interface/fluxx.h"
#include "array_macros/interface/fluxy.h"
#if NDIMS == 3
#include "array_macros/interface/fluxz.h"
#endif
#include "array_macros/interface/ifrcx.h"

#if NDIMS == 2
#define BEGIN \
  for(int cnt = 0, j = 1; j <= jsize; j++){ \
    for(int i = 2; i <= isize; i++, cnt++){
#define END \
    } \
  }
#else
#define BEGIN \
  for(int cnt = 0, k = 1; k <= ksize; k++){ \
    for(int j = 1; j <= jsize; j++){ \
      for(int i = 2; i <= isize; i++, cnt++){
#define END \
      } \
    } \
  }
#endif

static int advection_x(
    const domain_t * domain,
    const double * restrict ux,
    const double * restrict fluxx,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
#if NDIMS == 2
  BEGIN
    // ux is advected in x | 18
    const double hx_xm = HXXF(i-1);
    const double hx_x0 = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i-1);
    const double jd_x0 = JDXF(i  );
    const double jd_xp = JDXF(i+1);
    const double mux_xm = + 0.5 * jd_xm / hx_xm * FLUXX(i-1, j  )
                          + 0.5 * jd_x0 / hx_x0 * FLUXX(i  , j  );
    const double mux_xp = + 0.5 * jd_x0 / hx_x0 * FLUXX(i  , j  )
                          + 0.5 * jd_xp / hx_xp * FLUXX(i+1, j  );
    const double ux_xm = + 0.5 * UX(i-1, j  )
                         + 0.5 * UX(i  , j  );
    const double ux_xp = + 0.5 * UX(i  , j  )
                         + 0.5 * UX(i+1, j  );
    src[cnt] -= 1. / jd_x0 * (
        - mux_xm * ux_xm
        + mux_xp * ux_xp
    );
  END
#else
  BEGIN
    // ux is advected in x | 18
    const double hx_xm = HXXF(i-1);
    const double hx_x0 = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i-1);
    const double jd_x0 = JDXF(i  );
    const double jd_xp = JDXF(i+1);
    const double mux_xm = + 0.5 * jd_xm / hx_xm * FLUXX(i-1, j  , k  )
                          + 0.5 * jd_x0 / hx_x0 * FLUXX(i  , j  , k  );
    const double mux_xp = + 0.5 * jd_x0 / hx_x0 * FLUXX(i  , j  , k  )
                          + 0.5 * jd_xp / hx_xp * FLUXX(i+1, j  , k  );
    const double ux_xm = + 0.5 * UX(i-1, j  , k  )
                         + 0.5 * UX(i  , j  , k  );
    const double ux_xp = + 0.5 * UX(i  , j  , k  )
                         + 0.5 * UX(i+1, j  , k  );
    src[cnt] -= 1. / jd_x0 * (
        - mux_xm * ux_xm
        + mux_xp * ux_xp
    );
  END
#endif
  return 0;
}

static int advection_y(
    const domain_t * domain,
    const double * restrict ux,
    const double * restrict fluxy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  const double hy = domain->hy;
#if NDIMS == 2
  BEGIN
    // ux is advected in y | 15
    const double jd_xm = JDXC(i-1);
    const double jd_x0 = JDXF(i  );
    const double jd_xp = JDXC(i  );
    const double muy_ym = + 0.5 * jd_xm / hy * FLUXY(i-1, j  )
                          + 0.5 * jd_xp / hy * FLUXY(i  , j  );
    const double muy_yp = + 0.5 * jd_xm / hy * FLUXY(i-1, j+1)
                          + 0.5 * jd_xp / hy * FLUXY(i  , j+1);
    const double ux_ym = + 0.5 * UX(i  , j-1)
                         + 0.5 * UX(i  , j  );
    const double ux_yp = + 0.5 * UX(i  , j  )
                         + 0.5 * UX(i  , j+1);
    src[cnt] -= 1. / jd_x0 * (
        - muy_ym * ux_ym
        + muy_yp * ux_yp
    );
  END
#else
  BEGIN
    // ux is advected in y | 15
    const double jd_xm = JDXC(i-1);
    const double jd_x0 = JDXF(i  );
    const double jd_xp = JDXC(i  );
    const double muy_ym = + 0.5 * jd_xm / hy * FLUXY(i-1, j  , k  )
                          + 0.5 * jd_xp / hy * FLUXY(i  , j  , k  );
    const double muy_yp = + 0.5 * jd_xm / hy * FLUXY(i-1, j+1, k  )
                          + 0.5 * jd_xp / hy * FLUXY(i  , j+1, k  );
    const double ux_ym = + 0.5 * UX(i  , j-1, k  )
                         + 0.5 * UX(i  , j  , k  );
    const double ux_yp = + 0.5 * UX(i  , j  , k  )
                         + 0.5 * UX(i  , j+1, k  );
    src[cnt] -= 1. / jd_x0 * (
        - muy_ym * ux_ym
        + muy_yp * ux_yp
    );
  END
#endif
  return 0;
}

#if NDIMS == 3
static int advection_z(
    const domain_t * domain,
    const double * restrict ux,
    const double * restrict fluxz,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  const double hz = domain->hz;
  BEGIN
    // ux is advected in z | 15
    const double jd_xm = JDXC(i-1);
    const double jd_x0 = JDXF(i  );
    const double jd_xp = JDXC(i  );
    const double muz_zm = + 0.5 * jd_xm / hz * FLUXZ(i-1, j  , k  )
                          + 0.5 * jd_xp / hz * FLUXZ(i  , j  , k  );
    const double muz_zp = + 0.5 * jd_xm / hz * FLUXZ(i-1, j  , k+1)
                          + 0.5 * jd_xp / hz * FLUXZ(i  , j  , k+1);
    const double ux_zm = + 0.5 * UX(i  , j  , k-1)
                         + 0.5 * UX(i  , j  , k  );
    const double ux_zp = + 0.5 * UX(i  , j  , k  )
                         + 0.5 * UX(i  , j  , k+1);
    src[cnt] -= 1. / jd_x0 * (
        - muz_zm * ux_zm
        + muz_zp * ux_zp
    );
  END
  return 0;
}
#endif

static int diffusion_x(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict txx,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxc = domain->hxxc;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
#if NDIMS == 2
  BEGIN
    // ux is diffused in x | 11
    const double hx_xm = HXXC(i-1);
    const double hx_xp = HXXC(i  );
    const double jd_xm = JDXC(i-1);
    const double jd_x0 = JDXF(i  );
    const double jd_xp = JDXC(i  );
    const double txx_xm = TXX(i-1, j  );
    const double txx_xp = TXX(i  , j  );
    src[cnt] += diffusivity / jd_x0 * (
        - jd_xm / hx_xm * txx_xm
        + jd_xp / hx_xp * txx_xp
    );
  END
#else
  BEGIN
    // ux is diffused in x | 11
    const double hx_xm = HXXC(i-1);
    const double hx_xp = HXXC(i  );
    const double jd_xm = JDXC(i-1);
    const double jd_x0 = JDXF(i  );
    const double jd_xp = JDXC(i  );
    const double txx_xm = TXX(i-1, j  , k  );
    const double txx_xp = TXX(i  , j  , k  );
    src[cnt] += diffusivity / jd_x0 * (
        - jd_xm / hx_xm * txx_xm
        + jd_xp / hx_xp * txx_xp
    );
  END
#endif
  return 0;
}

static int diffusion_y(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict txy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double hy = domain->hy;
  const double * restrict jdxf = domain->jdxf;
#if NDIMS == 2
  BEGIN
    // ux is diffused in y | 7
    const double jd = JDXF(i  );
    const double txy_ym = TXY(i  , j  );
    const double txy_yp = TXY(i  , j+1);
    src[cnt] += diffusivity / jd * (
        - jd / hy * txy_ym
        + jd / hy * txy_yp
    );
  END
#else
  BEGIN
    // ux is diffused in y | 7
    const double jd = JDXF(i  );
    const double txy_ym = TXY(i  , j  , k  );
    const double txy_yp = TXY(i  , j+1, k  );
    src[cnt] += diffusivity / jd * (
        - jd / hy * txy_ym
        + jd / hy * txy_yp
    );
  END
#endif
  return 0;
}

#if NDIMS == 3
static int diffusion_z(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict txz,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict jdxf = domain->jdxf;
  BEGIN
    // ux is diffused in z | 7
    const double jd = JDXF(i  );
    const double txz_zm = TXZ(i  , j  , k  );
    const double txz_zp = TXZ(i  , j  , k+1);
    src[cnt] += diffusivity / jd * (
        - jd / hz * txz_zm
        + jd / hz * txz_zp
    );
  END
  return 0;
}
#endif

static int pressure(
    const domain_t * domain,
    const double * restrict p,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
#if NDIMS == 2
  BEGIN
    // pressure-gradient contribution | 4
    src[cnt] -= 1. / HXXF(i  ) * (
        - P(i-1, j  )
        + P(i  , j  )
    );
  END
#else
  BEGIN
    // pressure-gradient contribution | 4
    src[cnt] -= 1. / HXXF(i  ) * (
        - P(i-1, j  , k  )
        + P(i  , j  , k  )
    );
  END
#endif
  return 0;
}

static int surface(
    const domain_t * domain,
    const double * restrict ifrcx,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
#if NDIMS == 2
  BEGIN
    src[cnt] += IFRCX(i, j);
  END
#else
  BEGIN
    src[cnt] += IFRCX(i, j, k);
  END
#endif
  return 0;
}

static int gravity(
    const domain_t * domain,
    const double g,
    const double * den,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
#if NDIMS == 2
  BEGIN
    src[cnt] += g * (
        + 0.5 * DEN(i-1, j  )
        + 0.5 * DEN(i  , j  )
    );
  END
#else
  BEGIN
    src[cnt] += g * (
        + 0.5 * DEN(i-1, j  , k  )
        + 0.5 * DEN(i  , j  , k  )
    );
  END
#endif
  return 0;
}

/**
 * @brief comute right-hand-side of Runge-Kutta scheme of ux
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in,out] fluid  : n-step flow field (in), RK source terms (inout)
 * @return               : error code
 */
int compute_rhs_ux(
    const domain_t * domain,
    fluid_t * fluid,
    const interface_t * interface
){
  const double * restrict  ux = fluid-> ux.data;
  const double * restrict   p = fluid->  p.data;
  const double * restrict txx = fluid->txx.data;
  const double * restrict txy = fluid->txy.data;
#if NDIMS == 3
  const double * restrict txz = fluid->txz.data;
#endif
  const double * restrict fluxx = interface->fluxx.data;
  const double * restrict fluxy = interface->fluxy.data;
#if NDIMS == 3
  const double * restrict fluxz = interface->fluxz.data;
#endif
  double * restrict srca = fluid->srcux[rk_a].data;
  double * restrict srcg = fluid->srcux[rk_g].data;
  const double diffusivity = 1. / fluid->Re;
  const double acceleration = -1. / fluid->Fr / fluid->Fr;
  // advective contributions
  advection_x(domain, ux, fluxx, srca);
  advection_y(domain, ux, fluxy, srca);
#if NDIMS == 3
  advection_z(domain, ux, fluxz, srca);
#endif
  // diffusive contributions
  diffusion_x(domain, diffusivity, txx, srca);
  diffusion_y(domain, diffusivity, txy, srca);
#if NDIMS == 3
  diffusion_z(domain, diffusivity, txz, srca);
#endif
  // pressure-gradient contribution
  pressure(domain, p, srcg);
  surface(domain, interface->ifrcx.data, srca);
  gravity(domain, acceleration, fluid->den[0].data, srca);
  return 0;
}

/**
 * @brief predict ux
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : Runge-Kutta source terms (in), velocity (out)
 * @return               : error code
 */
int predict_ux(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double coef_a = rkcoefs[rkstep][rk_a];
  const double coef_b = rkcoefs[rkstep][rk_b];
  const double coef_g = rkcoefs[rkstep][rk_g];
  const double * restrict srcuxa = fluid->srcux[rk_a].data;
  const double * restrict srcuxb = fluid->srcux[rk_b].data;
  const double * restrict srcuxg = fluid->srcux[rk_g].data;
  double * restrict ux = fluid->ux.data;
  {
    const double * restrict den = fluid->den[0].data;
#if NDIMS == 2
    BEGIN
      const double lden = + 0.5 * DEN(i-1, j  )
                          + 0.5 * DEN(i  , j  );
      UX(i, j) =
        + lden * UX(i, j)
        + coef_a * dt * srcuxa[cnt]
        + coef_b * dt * srcuxb[cnt]
        + coef_g * dt * srcuxg[cnt];
    END
#else
    BEGIN
      const double lden = + 0.5 * DEN(i-1, j  , k  )
                          + 0.5 * DEN(i  , j  , k  );
      UX(i, j, k) =
        + lden * UX(i, j, k)
        + coef_a * dt * srcuxa[cnt]
        + coef_b * dt * srcuxb[cnt]
        + coef_g * dt * srcuxg[cnt];
    END
#endif
  }
  {
    const double * restrict den = fluid->den[1].data;
#if NDIMS == 2
    BEGIN
      const double lden = + 0.5 * DEN(i-1, j  )
                          + 0.5 * DEN(i  , j  );
      UX(i, j) /= lden;
    END
#else
    BEGIN
      const double lden = + 0.5 * DEN(i-1, j  , k  )
                          + 0.5 * DEN(i  , j  , k  );
      UX(i, j, k) /= lden;
    END
#endif
  }
  fluid_update_boundaries_ux(domain, &fluid->ux);
  return 0;
}

