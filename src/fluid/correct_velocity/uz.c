#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/psi.h"
#include "array_macros/fluid/den.h"

/**
 * @brief correct uz using scalar potential psi
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     dt_old : previous time step size
 * @param[in]     dt_new : current  time step size
 * @param[in,out] fluid  : scalar potential psi (in), ux (out)
 * @return               : error code
 */
int fluid_correct_velocity_uz(
    const domain_t * domain,
    const double dt_old,
    const double dt_new,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  double * restrict uz = fluid->uz.data;
  const double refden = fluid->refden;
  {
    const double * restrict psi = fluid->psi[1].data;
    // new scalar potential contribution
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          const double psi_zm = PSI(i  , j  , k-1);
          const double psi_zp = PSI(i  , j  , k  );
          UZ(i, j, k) -= dt_new / refden / hz * (
              - psi_zm
              + psi_zp
          );
        }
      }
    }
  }
  {
    const double * restrict psi = fluid->psi[0].data;
    const double * restrict den = fluid->den[1].data;
    const double coef = -1. / dt_new * dt_old;
    // old scalar potential contribution
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          const double psi_zm = PSI(i  , j  , k-1);
          const double psi_zp = PSI(i  , j  , k  );
          const double den_z0 = + 0.5 * DEN(i  , j  , k-1)
                                + 0.5 * DEN(i  , j  , k  );
          UZ(i, j, k) += dt_new * coef * (1. / den_z0 - 1. / refden) / hz * (
              - psi_zm
              + psi_zp
          );
        }
      }
    }
  }
  // update boundary and halo cells
  fluid_update_boundaries_uz(domain, &fluid->uz);
  return 0;
}
