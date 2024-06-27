#include <string.h>
#include "runge_kutta.h"
#include "array.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "interface.h"
#include "internal.h"

static int reset_srcs(
    const size_t rkstep,
    array_t * restrict srca,
    array_t * restrict srcb,
    array_t * restrict srcg
){
  // stash previous RK source term,
  //   which is achieved by swapping
  //   the pointers to "data"
  // NOTE: since "beta" is 0 when 0 == rkstep,
  //   this exchange is not needed
  if(0 != rkstep){
    double * tmp = srca->data;
    srca->data = srcb->data;
    srcb->data = tmp;
  }
  // zero-clear current RK source terms (exp/imp)
  memset(srca->data, 0, srca->datasize);
  memset(srcg->data, 0, srcg->datasize);
  return 0;
}

static int compute_rhs(
    const domain_t * domain,
    fluid_t * fluid,
    const interface_t * interface
){
  compute_rhs_ux(domain, fluid, interface);
  compute_rhs_uy(domain, fluid, interface);
  compute_rhs_uz(domain, fluid, interface);
  return 0;
}

static int predict(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
){
  predict_ux(domain, rkstep, dt, fluid);
  predict_uy(domain, rkstep, dt, fluid);
  predict_uz(domain, rkstep, dt, fluid);
  return 0;
}

/**
 * @brief predict the new velocity field and update the temperature field
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : flow field (in), RK source terms (in,out)
 * @return               : error code
 */
int fluid_predict_field(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid,
    const interface_t * interface
){
  // copy previous k-step source term and reset
  reset_srcs(rkstep, fluid->srcux + rk_a, fluid->srcux + rk_b, fluid->srcux + rk_g);
  reset_srcs(rkstep, fluid->srcuy + rk_a, fluid->srcuy + rk_b, fluid->srcuy + rk_g);
  reset_srcs(rkstep, fluid->srcuz + rk_a, fluid->srcuz + rk_b, fluid->srcuz + rk_g);
  // compute shear-stress tensor
  compute_txx(domain, fluid);
  compute_txy(domain, fluid);
  compute_txz(domain, fluid);
  compute_tyy(domain, fluid);
  compute_tyz(domain, fluid);
  compute_tzz(domain, fluid);
  // compute right-hand-side terms of the Runge-Kutta scheme
  compute_rhs(domain, fluid, interface);
  // update fields, which are still the prediction for the velocity,
  //   whereas the temperature is already updated to a new value
  predict(domain, rkstep, dt, fluid);
  return 0;
}

