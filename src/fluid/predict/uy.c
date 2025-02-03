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
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/p.h"
#include "array_macros/fluid/den.h"
#include "array_macros/fluid/txy.h"
#include "array_macros/fluid/tyy.h"
#include "array_macros/fluid/tyz.h"
#include "array_macros/interface/fluxx.h"
#include "array_macros/interface/fluxy.h"
#include "array_macros/interface/ifrcy.h"

#define BEGIN \
  for(int cnt = 0, j = 1; j <= jsize; j++){ \
    for(int i = 1; i <= isize; i++, cnt++){
#define END \
    } \
  }

static int advection_x(
    const domain_t * domain,
    const double * restrict uy,
    const double * restrict fluxx,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    // uy is advected in x
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
    const double mux_xm = + 0.5 * jd_xm / hx_xm * FLUXX(i  , j-1)
                          + 0.5 * jd_xm / hx_xm * FLUXX(i  , j  );
    const double mux_xp = + 0.5 * jd_xp / hx_xp * FLUXX(i+1, j-1)
                          + 0.5 * jd_xp / hx_xp * FLUXX(i+1, j  );
    const double uy_xm = + 0.5 * UY(i-1, j  )
                         + 0.5 * UY(i  , j  );
    const double uy_xp = + 0.5 * UY(i  , j  )
                         + 0.5 * UY(i+1, j  );
    src[cnt] -= 1. / jd_x0 * (
        - mux_xm * uy_xm
        + mux_xp * uy_xp
    );
  END
  return 0;
}

static int advection_y(
    const domain_t * domain,
    const double * restrict uy,
    const double * restrict fluxy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double hy = domain->hy;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    // uy is advected in y
    const double jd = JDXC(i  );
    const double muy_ym = + 0.5 * jd / hy * FLUXY(i  , j-1)
                          + 0.5 * jd / hy * FLUXY(i  , j  );
    const double muy_yp = + 0.5 * jd / hy * FLUXY(i  , j  )
                          + 0.5 * jd / hy * FLUXY(i  , j+1);
    const double uy_ym = + 0.5 * UY(i  , j-1)
                         + 0.5 * UY(i  , j  );
    const double uy_yp = + 0.5 * UY(i  , j  )
                         + 0.5 * UY(i  , j+1);
    src[cnt] -= 1. / jd * (
        - muy_ym * uy_ym
        + muy_yp * uy_yp
    );
  END
  return 0;
}

static int diffusion_x(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict txy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    // uy is diffused in x
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
    const double txy_xm = TXY(i  , j  );
    const double txy_xp = TXY(i+1, j  );
    src[cnt] += diffusivity / jd_x0 * (
        - jd_xm / hx_xm * txy_xm
        + jd_xp / hx_xp * txy_xp
    );
  END
  return 0;
}

static int diffusion_y(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict tyy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double hy = domain->hy;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    // uy is diffused in y
    const double jd = JDXC(i  );
    const double tyy_ym = TYY(i  , j-1);
    const double tyy_yp = TYY(i  , j  );
    src[cnt] += diffusivity / jd * (
        - jd / hy * tyy_ym
        + jd / hy * tyy_yp
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
  const double hy = domain->hy;
  BEGIN
    // pressure-gradient contribution
    src[cnt] -= 1. / hy * (
        - P(i  , j-1)
        + P(i  , j  )
    );
  END
  return 0;
}

static int surface(
    const domain_t * domain,
    const double * restrict ifrcy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  BEGIN
    src[cnt] += IFRCY(i, j);
  END
  return 0;
}

/**
 * @brief comute right-hand-side of Runge-Kutta scheme of uy
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in,out] fluid  : n-step flow field (in), RK source terms (inout)
 * @return               : error code
 */
int compute_rhs_uy(
    const domain_t * domain,
    fluid_t * fluid,
    const interface_t * interface
){
  const double * restrict  uy = fluid-> uy.data;
  const double * restrict   p = fluid->  p.data;
  const double * restrict txy = fluid->txy.data;
  const double * restrict tyy = fluid->tyy.data;
  const double * restrict fluxx = interface->fluxx.data;
  const double * restrict fluxy = interface->fluxy.data;
  double * restrict srca = fluid->srcuy[rk_a].data;
  double * restrict srcg = fluid->srcuy[rk_g].data;
  const double diffusivity = 1. / fluid->Re;
  // advective contributions
  advection_x(domain, uy, fluxx, srca);
  advection_y(domain, uy, fluxy, srca);
  // diffusive contributions
  diffusion_x(domain, diffusivity, txy, srca);
  diffusion_y(domain, diffusivity, tyy, srca);
  // pressure-gradient contribution
  pressure(domain, p, srcg);
  surface(domain, interface->ifrcy.data, srca);
  return 0;
}

/**
 * @brief predict uy
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : Runge-Kutta source terms (in), velocity (out)
 * @return               : error code
 */
int predict_uy(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double coef_a = rkcoefs[rkstep][rk_a];
  const double coef_b = rkcoefs[rkstep][rk_b];
  const double coef_g = rkcoefs[rkstep][rk_g];
  const double * restrict srcuya = fluid->srcuy[rk_a].data;
  const double * restrict srcuyb = fluid->srcuy[rk_b].data;
  const double * restrict srcuyg = fluid->srcuy[rk_g].data;
  double * restrict uy = fluid->uy.data;
  {
    const double * restrict den = fluid->den[0].data;
    BEGIN
      const double lden = + 0.5 * DEN(i  , j-1)
                          + 0.5 * DEN(i  , j  );
      UY(i, j) =
        + lden * UY(i, j)
        + coef_a * dt * srcuya[cnt]
        + coef_b * dt * srcuyb[cnt]
        + coef_g * dt * srcuyg[cnt];
    END
  }
  {
    const double * restrict den = fluid->den[1].data;
    BEGIN
      const double lden = + 0.5 * DEN(i  , j-1)
                          + 0.5 * DEN(i  , j  );
      UY(i, j) /= lden;
    END
  }
  fluid_update_boundaries_uy(domain, &fluid->uy);
  return 0;
}

