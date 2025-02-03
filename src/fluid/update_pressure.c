#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "array_macros/fluid/p.h"
#include "array_macros/fluid/psi.h"

static inline int add_explicit(
    const domain_t * domain,
    const fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict psi = fluid->psi[1].data;
  double * restrict p = fluid->p.data;
  // explicit contribution
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      P(i, j) += PSI(i, j);
    }
  }
  return 0;
}

static int update_scalar_potential(
    array_t * psi0,
    array_t * psi1
){
  // update psi
  double * tmp = psi0->data;
  psi0->data = psi1->data;
  psi1->data = tmp;
  return 0;
}

/**
 * @brief update pressure using scalar potential psi
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : scalar potential (in), pressure (out)
 * @return               : error code
 */
int fluid_update_pressure(
    const domain_t * domain,
    fluid_t * fluid
){
  // explicit contribution, always present
  add_explicit(domain, fluid);
  // impose boundary conditions and communicate halo cells
  fluid_update_boundaries_p(domain, &fluid->p);
  update_scalar_potential(&fluid->psi[0], &fluid->psi[1]);
  return 0;
}

