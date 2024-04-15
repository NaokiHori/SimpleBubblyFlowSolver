#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/psi.h"
#include "array_macros/fluid/den.h"

/**
 * @brief correct ux using scalar potential psi
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     dt_old : previous time step size
 * @param[in]     dt_new : current  time step size
 * @param[in,out] fluid  : scalar potential psi (in), ux (out)
 * @return               : error code
 */
int fluid_correct_velocity_ux(
    const domain_t * domain,
    const double dt_old,
    const double dt_new,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
  double * restrict ux = fluid->ux.data;
  const double refden = fluid->refden;
  {
    const double * restrict psi = fluid->psi[1].data;
#if NDIMS == 2
    // new scalar potential contribution | 11
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        const double hx = HXXF(i  );
        const double psi_xm = PSI(i-1, j  );
        const double psi_xp = PSI(i  , j  );
        UX(i, j) -= dt_new / refden / hx * (
            - psi_xm
            + psi_xp
        );
      }
    }
#else
    // new scalar potential contribution | 13
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 2; i <= isize; i++){
          const double hx = HXXF(i  );
          const double psi_xm = PSI(i-1, j  , k  );
          const double psi_xp = PSI(i  , j  , k  );
          UX(i, j, k) -= dt_new / refden / hx * (
              - psi_xm
              + psi_xp
          );
        }
      }
    }
#endif
  }
  {
    const double * restrict psi = fluid->psi[0].data;
    const double * restrict den = fluid->den[1].data;
    const double coef = -1. / dt_new * dt_old;
#if NDIMS == 2
    // old scalar potential contribution | 13
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        const double hx = HXXF(i  );
        const double psi_xm = PSI(i-1, j  );
        const double psi_xp = PSI(i  , j  );
        const double den_x0 = + 0.5 * DEN(i-1, j  )
                              + 0.5 * DEN(i  , j  );
        UX(i, j) += dt_new * coef * (1. / den_x0 - 1. / refden) / hx * (
            - psi_xm
            + psi_xp
        );
      }
    }
#else
    // old scalar potential contribution | 15
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 2; i <= isize; i++){
          const double hx = HXXF(i  );
          const double psi_xm = PSI(i-1, j  , k  );
          const double psi_xp = PSI(i  , j  , k  );
          const double den_x0 = + 0.5 * DEN(i-1, j  , k  )
                                + 0.5 * DEN(i  , j  , k  );
          UX(i, j, k) += dt_new * coef * (1. / den_x0 - 1. / refden) / hx * (
              - psi_xm
              + psi_xp
          );
        }
      }
    }
#endif
  }
  // update boundary and halo cells
  fluid_update_boundaries_ux(domain, &fluid->ux);
  return 0;
}

