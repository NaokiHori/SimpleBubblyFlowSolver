#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"

/**
 * @brief correct non-solenoidal velocity using scalar potential psi
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     dt_old : previous time step size
 * @param[in]     dt_new : current  time step size
 * @param[in,out] fluid  : scalar potential (in), velocity (out)
 * @return               : error code
 */
int fluid_correct_velocity(
    const domain_t * domain,
    const double dt_old,
    const double dt_new,
    fluid_t * fluid
){
  // compute prefactor gamma dt
  fluid_correct_velocity_ux(domain, dt_old, dt_new, fluid);
  fluid_correct_velocity_uy(domain, dt_old, dt_new, fluid);
#if NDIMS == 3
  fluid_correct_velocity_uz(domain, dt_old, dt_new, fluid);
#endif
  return 0;
}

