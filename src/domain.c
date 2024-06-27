#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include "sdecomp.h"
#include "memory.h"
#include "domain.h"
#include "fileio.h"
#include "array_macros/domain/xf.h"
#include "array_macros/domain/xc.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/hxxc.h"

/**
 * @brief load members in domain_t
 * @param[in]  dirname : name of directory from which data is loaded
 * @param[out] domain  : global domain sizes and resolutions
 * @return             : error code
 */
static int domain_load(
    const char dirname[],
    domain_t * domain
){
  size_t * glsizes = domain->glsizes;
  double * restrict lengths = domain->lengths;
  double * restrict * restrict xf = &domain->xf;
  double * restrict * restrict xc = &domain->xc;
  if(0 != fileio.r_serial(dirname, "glsizes", 1, (size_t [1]){NDIMS}, fileio.npy_size_t, sizeof(size_t), glsizes)){
    return 1;
  }
  if(0 != fileio.r_serial(dirname, "lengths", 1, (size_t [1]){NDIMS}, fileio.npy_double, sizeof(double), lengths)){
    return 1;
  }
  *xf = memory_calloc(glsizes[0] + 1, sizeof(double));
  if(0 != fileio.r_serial(dirname, "xf", 1, (size_t [1]){glsizes[0] + 1}, fileio.npy_double, sizeof(double), *xf)){
    return 1;
  }
  *xc = memory_calloc(glsizes[0] + 2, sizeof(double));
  if(0 != fileio.r_serial(dirname, "xc", 1, (size_t [1]){glsizes[0] + 2}, fileio.npy_double, sizeof(double), *xc)){
    return 1;
  }
  return 0;
}

/**
 * @brief save members in domain_t
 * @param[in] dirname : name of directory to which data is saved
 * @param[in] domain  : global domain sizes and resolutions
 * @return            : error code
 */
int domain_save(
    const char dirname[],
    const domain_t * domain
){
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  // since this is a serial operation,
  //   other processes are not involved
  if(root != myrank){
    return 0;
  }
  const size_t * glsizes = domain->glsizes;
  const double * lengths = domain->lengths;
  const double * xf      = domain->xf;
  const double * xc      = domain->xc;
  fileio.w_serial(dirname, "glsizes", 1, (size_t [1]){NDIMS}, fileio.npy_size_t, sizeof(size_t), glsizes);
  fileio.w_serial(dirname, "lengths", 1, (size_t [1]){NDIMS}, fileio.npy_double, sizeof(double), lengths);
  fileio.w_serial(dirname, "xf", 1, (size_t [1]){glsizes[0] + 1}, fileio.npy_double, sizeof(double), xf);
  fileio.w_serial(dirname, "xc", 1, (size_t [1]){glsizes[0] + 2}, fileio.npy_double, sizeof(double), xc);
  return 0;
}

static int check_x_grid_is_uniform(
    const domain_t * domain,
    bool * x_grid_is_uniform
){
  static bool is_uniform = false;
  const size_t isize = domain->glsizes[0];
  const double * hxxc = domain->hxxc;
  double extrema[2] = {+1. * DBL_MAX, -1. * DBL_MAX};
  for(size_t i = 0; i < isize; i++){
    extrema[0] = fmin(extrema[0], hxxc[i]);
    extrema[1] = fmax(extrema[1], hxxc[i]);
  }
  if(fabs(extrema[1] - extrema[0]) < 1.e-15){
    is_uniform = true;
  }else{
    is_uniform = false;
  }
  *x_grid_is_uniform = is_uniform;
  return 0;
}

// scale factor in x, defined at x cell faces
static double * allocate_and_init_hxxf(
    const int isize,
    const double * xc
){
  double * hxxf = memory_calloc(isize + 1, sizeof(double));
  for(int i = 1; i <= isize + 1; i++){
    HXXF(i  ) = XC(i  ) - XC(i-1);
  }
  return hxxf;
}

// scale factor in x, defined at x cell centers
static double * allocate_and_init_hxxc(
    const int isize,
    const double * xf
){
  double * hxxc = memory_calloc(isize, sizeof(double));
  for(int i = 1; i <= isize; i++){
    HXXC(i  ) = XF(i+1) - XF(i  );
  }
  return hxxc;
}

// Jacobian determinant, defined at x cell faces
static double * allocate_and_init_jdxf(
    const int isize,
    const double * hxxf,
    const double hy
){
  double * jdxf = memory_calloc(isize + 1, sizeof(double));
  for(int i = 1; i <= isize + 1; i++){
    JDXF(i  ) = HXXF(i  ) * hy;
  }
  return jdxf;
}

// Jacobian determinant, defined at x cell centers
static double * allocate_and_init_jdxc(
    const int isize,
    const double * hxxc,
    const double hy
){
  double * jdxc = memory_calloc(isize, sizeof(double));
  for(int i = 1; i <= isize; i++){
    JDXC(i  ) = HXXC(i  ) * hy;
  }
  return jdxc;
}

static void report(
    const domain_t * domain
){
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(root == myrank){
    printf("DOMAIN\n");
    for(sdecomp_dir_t dim = 0; dim < NDIMS; dim++){
      printf("\tglsizes[%u]: %zu\n", dim, domain->glsizes[dim]);
    }
    for(sdecomp_dir_t dim = 0; dim < NDIMS; dim++){
      printf("\tlengths[%u]: % .7e\n", dim, domain->lengths[dim]);
    }
    fflush(stdout);
  }
}

/**
 * @brief constructor of the structure
 * @param[in]  dirname_ic : name of directory in which initial conditions are stored
 * @param[out] domain     : structure being allocated and initalised
 * @return                : (success) 0
 *                          (failure) non-zero value
 */
int domain_init(
    const char dirname_ic[],
    domain_t * domain
){
  sdecomp_info_t ** info    = &domain->info;
  size_t * restrict glsizes =  domain->glsizes;
  size_t * restrict mysizes =  domain->mysizes;
  size_t * restrict offsets =  domain->offsets;
  double * restrict lengths =  domain->lengths;
  double * restrict * xf    = &domain->xf;
  double * restrict * xc    = &domain->xc;
  double * restrict * jdxf  = &domain->jdxf;
  double * restrict * jdxc  = &domain->jdxc;
  double * restrict * hxxf   = &domain->hxxf;
  double * restrict * hxxc   = &domain->hxxc;
  double * restrict   hy    = &domain->hy;
  // load spatial information
  if(0 != domain_load(dirname_ic, domain)){
    return 1;
  }
  // scale factors
  *hxxf = allocate_and_init_hxxf(glsizes[0], *xc);
  *hxxc = allocate_and_init_hxxc(glsizes[0], *xf);
  *hy = lengths[1] / glsizes[1];
  // Jacobian determinants at cell faces and centers
  *jdxf = allocate_and_init_jdxf(glsizes[0], *hxxf, *hy);
  *jdxc = allocate_and_init_jdxc(glsizes[0], *hxxc, *hy);
  // initialise sdecomp to distribute the domain
  if(0 != sdecomp.construct(
        MPI_COMM_WORLD,
        NDIMS,
        (size_t [NDIMS]){0, 0},
        (bool [NDIMS]){false, true},
        info
  )) return 1;
  // local array sizes and offsets
  for(size_t dim = 0; dim < NDIMS; dim++){
    sdecomp.get_pencil_mysize(*info, SDECOMP_X1PENCIL, dim, glsizes[dim], mysizes + dim);
    sdecomp.get_pencil_offset(*info, SDECOMP_X1PENCIL, dim, glsizes[dim], offsets + dim);
  }
  // in this project, only uniform grids are allowed
  bool is_uniform = false;
  check_x_grid_is_uniform(domain, &is_uniform);
  if (!is_uniform) {
    fprintf(stderr, "x grid is not uniform, which is not allowed");
    return 1;
  }
  report(domain);
  return 0;
}

