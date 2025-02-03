#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/psi.h"
#include "array_macros/fluid/den.h"

/**
 * @brief correct uy using scalar potential psi
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     dt_old : previous time step size
 * @param[in]     dt_new : current  time step size
 * @param[in,out] fluid  : scalar potential psi (in), ux (out)
 * @return               : error code
 */
int fluid_correct_velocity_uy(
    const domain_t * domain,
    const double dt_old,
    const double dt_new,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double hy = domain->hy;
  double * restrict uy = fluid->uy.data;
  const double refden = fluid->refden;
  // new scalar potential
  {
    const double * restrict psi = fluid->psi[1].data;
    // new scalar potential contribution
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double psi_ym = PSI(i  , j-1);
        const double psi_yp = PSI(i  , j  );
        UY(i, j) -= dt_new / refden / hy * (
            - psi_ym
            + psi_yp
        );
      }
    }
  }
  {
    const double * restrict psi = fluid->psi[0].data;
    const double * restrict den = fluid->den[1].data;
    const double coef = -1. / dt_new * dt_old;
    // old scalar potential contribution
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double psi_ym = PSI(i  , j-1);
        const double psi_yp = PSI(i  , j  );
        const double den_y0 = + 0.5 * DEN(i  , j-1)
                              + 0.5 * DEN(i  , j  );
        UY(i, j) += dt_new * coef * (1. / den_y0 - 1. / refden) / hy * (
            - psi_ym
            + psi_yp
        );
      }
    }
  }
  // update boundary and halo cells
  fluid_update_boundaries_uy(domain, &fluid->uy);
  return 0;
}

