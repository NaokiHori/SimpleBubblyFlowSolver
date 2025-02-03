#include <math.h>
#include "memory.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "fileio.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/p.h"
#include "array_macros/fluid/den.h"
#include "array_macros/fluid/visux.h"
#include "array_macros/fluid/visuy.h"
#include "array_macros/fluid/visuz.h"
#include "array_macros/fluid/txx.h"
#include "array_macros/fluid/txy.h"
#include "array_macros/fluid/txz.h"
#include "array_macros/fluid/tyy.h"
#include "array_macros/fluid/tyz.h"
#include "array_macros/fluid/tzz.h"
#include "array_macros/fluid/psi.h"
#include "array_macros/fluid/srcux.h"
#include "array_macros/fluid/srcuy.h"
#include "array_macros/fluid/srcuz.h"

/**
 * @brief allocate members
 * @param[in]  domain : information about domain decomposition and size
 * @param[out] fluid  : structure storing flow fields and auxiliary buffers
 * @return            : error code
 */
static int allocate(
    const domain_t * domain,
    fluid_t * fluid
){
  // velocity
  if(0 != array.prepare(domain, UX_NADDS, sizeof(double), &fluid->ux )) return 1;
  if(0 != array.prepare(domain, UY_NADDS, sizeof(double), &fluid->uy )) return 1;
  if(0 != array.prepare(domain, UZ_NADDS, sizeof(double), &fluid->uz )) return 1;
  // pressure and scalar potentials
  if(0 != array.prepare(domain, P_NADDS,   sizeof(double), &fluid->p  )) return 1;
  for(size_t n = 0; n < 2; n++){
    if(0 != array.prepare(domain, PSI_NADDS, sizeof(double), &fluid->psi[n])) return 1;
  }
  // density
  for(size_t n = 0; n < 2; n++){
    if(0 != array.prepare(domain, DEN_NADDS, sizeof(double), &fluid->den[n])) return 1;
  }
  // viscosity
  if(0 != array.prepare(domain, VISUX_NADDS, sizeof(double), &fluid->visux)) return 1;
  if(0 != array.prepare(domain, VISUY_NADDS, sizeof(double), &fluid->visuy)) return 1;
  if(0 != array.prepare(domain, VISUZ_NADDS, sizeof(double), &fluid->visuz)) return 1;
  // stress tensor
  if(0 != array.prepare(domain, TXX_NADDS, sizeof(double), &fluid->txx)) return 1;
  if(0 != array.prepare(domain, TXY_NADDS, sizeof(double), &fluid->txy)) return 1;
  if(0 != array.prepare(domain, TXZ_NADDS, sizeof(double), &fluid->txz)) return 1;
  if(0 != array.prepare(domain, TYY_NADDS, sizeof(double), &fluid->tyy)) return 1;
  if(0 != array.prepare(domain, TYZ_NADDS, sizeof(double), &fluid->tyz)) return 1;
  if(0 != array.prepare(domain, TZZ_NADDS, sizeof(double), &fluid->tzz)) return 1;
  // Runge-Kutta source terms
  for(size_t n = 0; n < 3; n++){
    if(0 != array.prepare(domain, SRCUX_NADDS, sizeof(double), &fluid->srcux[n])) return 1;
    if(0 != array.prepare(domain, SRCUY_NADDS, sizeof(double), &fluid->srcuy[n])) return 1;
    if(0 != array.prepare(domain, SRCUZ_NADDS, sizeof(double), &fluid->srcuz[n])) return 1;
  }
  return 0;
}

static void report(
    const sdecomp_info_t * info,
    const fluid_t * fluid
){
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(info, &myrank);
  if(root == myrank){
    printf("FLUID\n");
    printf("\tRe: % .7e\n", fluid->Re);
    printf("\tFr: % .7e\n", fluid->Fr);
    printf("\tDensity ratio: % .7e\n", fluid->denr);
    printf("\tViscosity ratio: % .7e\n", fluid->visr);
    fflush(stdout);
  }
}

/**
 * @brief constructor of the structure
 * @param[in]  dirname_ic : name of directory in which initial flow fields are stored
 * @param[in]  domain     : information about domain decomposition and size
 * @param[out]            : structure being allocated and initalised
 * @return                : (success) 0
 *                          (failure) non-zero value
 */
int fluid_init(
    const char dirname_ic[],
    const domain_t * domain,
    fluid_t * fluid
){
  // allocate arrays
  if(0 != allocate(domain, fluid)) return 1;
  // load flow fields
  if(0 != array.load(domain, dirname_ic,  "ux", fileio.npy_double, &fluid->    ux)) return 1;
  if(0 != array.load(domain, dirname_ic,  "uy", fileio.npy_double, &fluid->    uy)) return 1;
  if(0 != array.load(domain, dirname_ic,  "uz", fileio.npy_double, &fluid->    uz)) return 1;
  if(0 != array.load(domain, dirname_ic,   "p", fileio.npy_double, &fluid->     p)) return 1;
  if(0 != array.load(domain, dirname_ic, "psi", fileio.npy_double, &fluid->psi[0])) return 1;
  // impose boundary conditions and communicate halo cells
  fluid_update_boundaries_ux (domain, &fluid->    ux);
  fluid_update_boundaries_uy (domain, &fluid->    uy);
  fluid_update_boundaries_uz (domain, &fluid->    uz);
  fluid_update_boundaries_p  (domain, &fluid->     p);
  fluid_update_boundaries_psi(domain, &fluid->psi[0]);
  // compute diffusivities
  if(0 != config.get_double("Re", &fluid->Re)) return 1;
  if(0 != config.get_double("Fr", &fluid->Fr)) return 1;
  if(0 != config.get_double("denr", &fluid->denr)) return 1;
  if(0 != config.get_double("visr", &fluid->visr)) return 1;
  fluid->refden = fmin(1., fluid->denr);
  report(domain->info, fluid);
  return 0;
}

