#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include "config.h"
#include "array.h"
#include "sdecomp.h"
#include "domain.h"
#include "fluid.h"
#include "interface.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/domain/hxxf.h"

static const double pi = 3.1415926535897932384626433832;

// overriden later using environment variables
static bool initialised = false;
static double coef_dt_adv = 0.;
static double coef_dt_dif = 0.;
static double coef_dt_int = 0.;

// max possible dt
static const double dt_max = 1.;

/**
 * @brief decide time step size restricted by the advective terms
 * @param[in]  domain : information about domain decomposition and size
 * @param[in]  fluid  : velocity
 * @param[out] dt     : time step size
 * @return            : error code
 */
static int decide_dt_adv(
    const domain_t * domain,
    const fluid_t * fluid,
    double * restrict dt
){
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double hy = domain->hy;
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
  // sufficiently small number to avoid zero division
  const double small = 1.e-8;
  *dt = dt_max;
  // compute grid-size over velocity in x
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      const double hx = HXXF(i  );
      const double vel = fabs(UX(i, j)) + small;
      *dt = fmin(*dt, hx / vel);
    }
  }
  // compute grid-size over velocity in y
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      const double vel = fabs(UY(i, j)) + small;
      *dt = fmin(*dt, hy / vel);
    }
  }
  // compute grid-size over velocity in z
  // unify result, multiply safety factor
  MPI_Allreduce(MPI_IN_PLACE, dt, 1, MPI_DOUBLE, MPI_MIN, comm_cart);
  *dt *= coef_dt_adv;
  return 0;
}

/**
 * @brief decide time step size restricted by the diffusive terms
 * @param[in]  domain      : grid size
 * @param[in]  diffusivity : fluid / temperature diffusivity
 * @param[out] dt          : time step size
 * @return                 : error code
 */
static int decide_dt_dif(
    const domain_t * domain,
    const double denr,
    const double visr,
    const double diffusivity,
    double * restrict dt
){
  const int isize = domain->mysizes[0];
  const double * restrict hxxf = domain->hxxf;
  const double hy = domain->hy;
  double grid_sizes[NDIMS] = {0.};
  // find minimum grid size in x direction
  grid_sizes[0] = DBL_MAX;
  for(int i = 2; i <= isize; i++){
    const double hx = HXXF(i  );
    grid_sizes[0] = fmin(grid_sizes[0], hx);
  }
  grid_sizes[1] = hy;
  // compute diffusive constraints
  for(int dim = 0; dim < NDIMS; dim++){
    dt[dim] = fmin(1., denr / visr) / diffusivity * 0.5 / NDIMS * pow(grid_sizes[dim], 2.);
  }
  for(int dim = 0; dim < NDIMS; dim++){
    dt[dim] *= coef_dt_dif;
  }
  return 0;
}

static int decide_dt_int(
    const domain_t * domain,
    const double tension,
    double * restrict dt
){
  const int isize = domain->mysizes[0];
  const double * restrict hxxf = domain->hxxf;
  const double hy = domain->hy;
  double grid_sizes[NDIMS] = {0.};
  // find minimum grid size in x direction
  grid_sizes[0] = DBL_MAX;
  for(int i = 2; i <= isize; i++){
    const double hx = HXXF(i  );
    grid_sizes[0] = fmin(grid_sizes[0], hx);
  }
  grid_sizes[1] = hy;
  // compute interfacial constraints
  *dt = dt_max;
  for(int dim = 0; dim < NDIMS; dim++){
    *dt = fmin(*dt, sqrt(1. / tension / 4. / pi * pow(grid_sizes[dim], 3.)));
  }
  *dt *= coef_dt_int;
  return 0;
}

/**
 * @brief decide time step size which can integrate the equations stably
 * @param[in]  domain : information about domain decomposition and size
 * @param[in]  fluid  : velocity and diffusivities
 * @param[out] dt     : time step size
 * @return            : (success) 0
 *                    : (failure) non-zero value
 */
int decide_dt(
    const domain_t * domain,
    const fluid_t * fluid,
    const interface_t * interface,
    double * dt
){
  // load safety factors
  if(!initialised){
    if(0 != config.get_double("coef_dt_adv", &coef_dt_adv)) return 1;
    if(0 != config.get_double("coef_dt_dif", &coef_dt_dif)) return 1;
    if(0 != config.get_double("coef_dt_int", &coef_dt_int)) return 1;
    initialised = true;
  }
  // compute advective and diffusive constraints
  double dt_adv[    1] = {0.};
  double dt_dif[NDIMS] = {0.};
  double dt_int[    1] = {0.};
  decide_dt_adv(domain, fluid, dt_adv);
  decide_dt_dif(domain, fluid->denr, fluid->visr, 1. / fluid->Re, dt_dif);
  decide_dt_int(domain, 1. / interface->We, dt_int);
  // choose smallest value as dt
  // advection
  *dt = dt_adv[0];
  // diffusion, momentum
  *dt = fmin(*dt, dt_dif[0]);
  *dt = fmin(*dt, dt_dif[1]);
  // surface tension
  *dt = fmin(*dt, dt_int[0]);
  return 0;
}

