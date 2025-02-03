#include <math.h>
#include <mpi.h>
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "interface.h"
#include "array_macros/fluid/visux.h"
#include "array_macros/fluid/visuy.h"
#include "array_macros/fluid/visuz.h"
#include "array_macros/interface/vof.h"

static inline double get(
    const double visr,
    const double min,
    const double max,
    const double vof
){
  double vis = 1. + (visr - 1.) * vof;
  vis = fmax(min, vis);
  vis = fmin(max, vis);
  return vis;
}

static int compute_x(
    const domain_t * domain,
    const array_t * restrict arr_vof,
    const double visr,
    const double min,
    const double max,
    array_t * restrict arr_visux
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict vof = arr_vof->data;
  double * restrict visux = arr_visux->data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      VISUX(1, j, k) = get(
          visr, min, max,
          VOF(0, j, k)
      );
      for(int i = 2; i <= isize; i++){
        VISUX(i, j, k) = get(
            visr, min, max,
            + 0.5 * VOF(i-1, j  , k  )
            + 0.5 * VOF(i  , j  , k  )
        );
      }
      VISUX(isize + 1, j, k) = get(
          visr, min, max,
          VOF(isize + 1, j, k)
      );
    }
  }
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, arr_visux)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, arr_visux)){
    return 1;
  }
  return 0;
}

static int compute_y(
    const domain_t * domain,
    const array_t * restrict arr_vof,
    const double visr,
    const double min,
    const double max,
    array_t * restrict arr_visuy
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict vof = arr_vof->data;
  double * restrict visuy = arr_visuy->data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize + 1; i++){
        VISUY(i, j, k) = get(
            visr, min, max,
            + 0.5 * VOF(i  , j-1, k  )
            + 0.5 * VOF(i  , j  , k  )
        );
      }
    }
  }
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, arr_visuy)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, arr_visuy)){
    return 1;
  }
  return 0;
}

static int compute_z(
    const domain_t * domain,
    const array_t * restrict arr_vof,
    const double visr,
    const double min,
    const double max,
    array_t * restrict arr_visuz
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict vof = arr_vof->data;
  double * restrict visuz = arr_visuz->data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize + 1; i++){
        VISUZ(i, j, k) = get(
            visr, min, max,
            + 0.5 * VOF(i  , j  , k-1)
            + 0.5 * VOF(i  , j  , k  )
        );
      }
    }
  }
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, arr_visuz)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, arr_visuz)){
    return 1;
  }
  return 0;
}

int fluid_compute_viscosity(
    const domain_t * domain,
    fluid_t * fluid,
    const interface_t * interface
){
  const double visr = fluid->visr;
  const double min = fmin(1., visr);
  const double max = fmax(1., visr);
  if(0 != compute_x(domain, &interface->vof, visr, min, max, &fluid->visux)){
    return 1;
  }
  if(0 != compute_y(domain, &interface->vof, visr, min, max, &fluid->visuy)){
    return 1;
  }
  if(0 != compute_z(domain, &interface->vof, visr, min, max, &fluid->visuz)){
    return 1;
  }
  return 0;
}

