#include "config.h"
#include "domain.h"
#include "interface.h"
#include "interface_solver.h"
#include "fileio.h"
#include "array_macros/interface/vof.h"
#include "array_macros/interface/ifrcx.h"
#include "array_macros/interface/ifrcy.h"
#include "array_macros/interface/ifrcz.h"
#include "array_macros/interface/dvof.h"
#include "array_macros/interface/normal.h"
#include "array_macros/interface/curv.h"
#include "array_macros/interface/fluxx.h"
#include "array_macros/interface/fluxy.h"
#include "array_macros/interface/fluxz.h"
#include "array_macros/interface/src.h"

/**
 * @brief allocate interface_t
 * @param[in]  domain    : information about domain decomposition and size
 * @param[out] interface : structure being allocated
 * @return               : error code
 */
static int allocate(
    const domain_t * domain,
    interface_t * interface
){
  if(0 != array.prepare(domain, VOF_NADDS, sizeof(double), &interface->vof)) return 1;
  if(0 != array.prepare(domain, IFRCX_NADDS, sizeof(double), &interface->ifrcx)) return 1;
  if(0 != array.prepare(domain, IFRCY_NADDS, sizeof(double), &interface->ifrcy)) return 1;
  if(0 != array.prepare(domain, IFRCZ_NADDS, sizeof(double), &interface->ifrcz)) return 1;
  if(0 != array.prepare(domain, DVOF_NADDS, sizeof(vector_t), &interface->dvof)) return 1;
  if(0 != array.prepare(domain, NORMAL_NADDS, sizeof(normal_t), &interface->normal)) return 1;
  if(0 != array.prepare(domain, CURV_NADDS, sizeof(double), &interface->curv)) return 1;
  if(0 != array.prepare(domain, FLUXX_NADDS, sizeof(double), &interface->fluxx)) return 1;
  if(0 != array.prepare(domain, FLUXY_NADDS, sizeof(double), &interface->fluxy)) return 1;
  if(0 != array.prepare(domain, FLUXZ_NADDS, sizeof(double), &interface->fluxz)) return 1;
  for(size_t n = 0; n < 2; n++){
    if(0 != array.prepare(domain, SRC_NADDS, sizeof(double), &interface->src[n])) return 1;
  }
  return 0;
}

static void report(
    const sdecomp_info_t * info,
    const interface_t * interface
){
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(info, &myrank);
  if(root == myrank){
    FILE * stream = stdout;
    fprintf(stream, "INTERFACE\n");
    fprintf(stream, "\tWe: % .7e\n", interface->We);
    fflush(stream);
  }
}

int interface_init(
    const char dirname_ic[],
    const domain_t * domain,
    interface_t * interface
){
  if(0 != allocate(domain, interface)) return 1;
  // load interface field
  if(0 != array.load(domain, dirname_ic, "vof", fileio.npy_double, &interface->vof)) return 1;
  // impose boundary conditions and communicate halo cells
  interface_update_boundaries_vof(domain, &interface->vof);
  // compute surface tension coefficient
  if(0 != config.get_double("We", &interface->We)) return 1;
  report(domain->info, interface);
  return 0;
}

