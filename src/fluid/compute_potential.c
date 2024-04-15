#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "sdecomp.h"
#include "memory.h"
#include "domain.h"
#include "tdm.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#if NDIMS == 3
#include "array_macros/fluid/uz.h"
#endif
#include "array_macros/fluid/den.h"
#include "array_macros/fluid/psi.h"

static const double g_pi = 3.14159265358979324;

// structure only used to solve Poisson equation
// NOTE: this is shared among normal and efficient solvers
//   thus some variables may not be used
typedef struct {
  bool is_initialised;
  void * restrict buf0;
  void * restrict buf1;
  fftw_plan fftw_plan_x[2];
#if NDIMS == 3
  fftw_plan fftw_plan_y[2];
  fftw_plan fftw_plan_z[2];
#endif
  size_t tdm_sizes[2];
  tdm_info_t * tdm_info;
  double * evals;
  sdecomp_transpose_plan_t * r_transposer_x1_to_y1;
  sdecomp_transpose_plan_t * r_transposer_y1_to_x1;
#if NDIMS == 3
  sdecomp_transpose_plan_t * c_transposer_y1_to_z1;
  sdecomp_transpose_plan_t * c_transposer_z1_to_y1;
#endif
} poisson_solver_t;

/* initialise Poisson solver */
// several pencils for different data types are treated
//   and thus this source is very complicated
// used prefixes are as follows:
//   r_: real    (double)       type
//   c_: complex (fftw_complex) type
//
//   gl_    : global array size (not pencils)
//   x1pncl_: each pencl (x1, y1, ...)

/* size of domain and pencils */
// NOTE: define globally to reduce the number of arguments of static functions
// global domain size in real space
static size_t r_gl_sizes[NDIMS] = {0};
#if NDIMS == 3
// global domain size in complex space
static size_t c_gl_sizes[NDIMS] = {0};
#endif
// local domain size (x1 pencil) in real space
static size_t r_x1pncl_sizes[NDIMS] = {0};
// local domain size (y1 pencil) in real space
static size_t r_y1pncl_sizes[NDIMS] = {0};
#if NDIMS == 3
// local domain size (y1 pencil) in complex space
static size_t c_y1pncl_sizes[NDIMS] = {0};
// local domain size (z1 pencil) in complex space
static size_t c_z1pncl_sizes[NDIMS] = {0};
#endif

static size_t prod(
    const size_t sizes[NDIMS]
){
  // compute the product of the given vector
  size_t nitems = 1;
  for(size_t dim = 0; dim < NDIMS; dim++){
    nitems *= sizes[dim];
  }
  return nitems;
}

static int report_failure(
    const char type[]
){
  // function to just dump error message and abort
  FILE * stream = stderr;
  fprintf(stream, "Poisson solver, initialisation failed: %s\n", type);
  fprintf(stream, "  FFTW:    A possible reason is you link Intel-MKL lib\n");
  fprintf(stream, "           Make sure you use FFTW3 directly,\n");
  fprintf(stream, "           NOT its crazy wrapper offered by MKL\n");
  fprintf(stream, "  SDECOMP: Check sdecomp.log and check arguments\n");
  fprintf(stream, "           If they are all correct, PLEASE CONTACT ME\n");
  fflush(stream);
  return 0;
}

static size_t max(
    const size_t val0,
    const size_t val1
){
  if(val0 > val1){
    return val0;
  }else{
    return val1;
  }
}

static int compute_pencil_sizes(
    const domain_t * domain
){
  // NOTE: those variables are defined globally at the top of this file
  //   to reduce the nhumber of arguments which functions take
  // global domain size in real space
  const sdecomp_info_t * info = domain->info;
  r_gl_sizes[0] = domain->glsizes[0];
  r_gl_sizes[1] = domain->glsizes[1];
#if NDIMS == 3
  r_gl_sizes[2] = domain->glsizes[2];
#endif
#if NDIMS == 3
  // global domain size in complex space
  // NOTE: Hermite symmetry in y
  c_gl_sizes[0] = domain->glsizes[0];
  c_gl_sizes[1] = domain->glsizes[1] / 2 + 1;
  c_gl_sizes[2] = domain->glsizes[2];
#endif
  // local domain sizes
  for(sdecomp_dir_t dim = 0; dim < NDIMS; dim++){
    if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_X1PENCIL, dim, r_gl_sizes[dim], r_x1pncl_sizes + dim)) return 1;
    if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, r_gl_sizes[dim], r_y1pncl_sizes + dim)) return 1;
#if NDIMS == 3
    if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, c_gl_sizes[dim], c_y1pncl_sizes + dim)) return 1;
    if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_Z1PENCIL, dim, c_gl_sizes[dim], c_z1pncl_sizes + dim)) return 1;
#endif
  }
  return 0;
}

static int allocate_buffers(
    poisson_solver_t * poisson_solver
){
  // although there are bunch of pencils involved,
  //   two buffers are enough to do the job,
  //   which are allocated here
  void * restrict * buf0 = &poisson_solver->buf0;
  void * restrict * buf1 = &poisson_solver->buf1;
  const size_t r_dsize = sizeof(double);
#if NDIMS == 3
  const size_t c_dsize = sizeof(fftw_complex);
#endif
  size_t buf0_bytes = 0;
  size_t buf1_bytes = 0;
#if NDIMS == 2
  // r_x1pncl -> FFT -> r_x1pncl -> rotate -> r_y1pncl
  // buffer0            buffer1               buffer0
  buf0_bytes = max(buf0_bytes, r_dsize * prod(r_x1pncl_sizes));
  buf0_bytes = max(buf0_bytes, r_dsize * prod(r_y1pncl_sizes));
  buf1_bytes = max(buf1_bytes, r_dsize * prod(r_x1pncl_sizes));
#else
  // r_x1pncl -> FFT -> r_x1pncl -> rotate -> r_y1pncl -> FFT -> c_y1pncl -> rotate -> c_z1pncl
  // buffer0            buffer1               buffer0            buffer1               buffer0
  buf0_bytes = max(buf0_bytes, r_dsize * prod(r_x1pncl_sizes));
  buf0_bytes = max(buf0_bytes, r_dsize * prod(r_y1pncl_sizes));
  buf0_bytes = max(buf0_bytes, c_dsize * prod(c_z1pncl_sizes));
  buf1_bytes = max(buf1_bytes, r_dsize * prod(r_x1pncl_sizes));
  buf1_bytes = max(buf1_bytes, c_dsize * prod(c_y1pncl_sizes));
#endif
  // allocate them using fftw_malloc to enforce them 16bit-aligned for SIMD
  *buf0 = fftw_malloc(buf0_bytes);
  if(NULL == *buf0){
    fprintf(stderr, "FATAL: fftw_malloc failed (requested %zu bytes)\n", buf0_bytes);
    fflush(stderr);
    return 1;
  }
  *buf1 = fftw_malloc(buf1_bytes);
  if(NULL == *buf1){
    fprintf(stderr, "FATAL: fftw_malloc failed (requested %zu bytes)\n", buf1_bytes);
    fflush(stderr);
    return 1;
  }
  return 0;
}

static int init_tri_diagonal_solver(
    const domain_t * domain,
    poisson_solver_t * poisson_solver
){
  // N x N tri-diagonal matrix,
  //   which are solved for M times
  // tdm_sizes[0] = N, tdm_sizes[1] = M
  // since lower- and upper-diagonal components are
  //   independent to y, z directions and time,
  //   we compute here and re-use them
  // center-diagonal components are, on the other hand,
  //   dependent on time and thus needs to compute everytime
  //   in the solver
  size_t * restrict tdm_sizes = poisson_solver->tdm_sizes;
  tdm_info_t ** tdm_info = &poisson_solver->tdm_info;
#if NDIMS == 2
  tdm_sizes[0] = r_y1pncl_sizes[1];
  tdm_sizes[1] = r_y1pncl_sizes[0];
  if(0 != tdm.construct(
    /* size of system */ tdm_sizes[0],
    /* number of rhs  */ 1,
    /* is periodic    */ true,
    /* is complex     */ false,
    /* output         */ tdm_info
  )) return 1;
  // initialise tri-diagonal matrix in y direction
  double * tdm_l = NULL;
  double * tdm_u = NULL;
  tdm.get_l(*tdm_info, &tdm_l);
  tdm.get_u(*tdm_info, &tdm_u);
  const double hy = domain->hy;
  for(size_t j = 0; j < tdm_sizes[0]; j++){
    tdm_l[j] = 1. / hy / hy;
    tdm_u[j] = 1. / hy / hy;
  }
#else
  tdm_sizes[0] = c_z1pncl_sizes[2];
  tdm_sizes[1] = c_z1pncl_sizes[0] * c_z1pncl_sizes[1];
  if(0 != tdm.construct(
    /* size of system */ tdm_sizes[0],
    /* number of rhs  */ 1,
    /* is periodic    */ true,
    /* is complex     */ true,
    /* output         */ tdm_info
  )) return 1;
  // initialise tri-diagonal matrix in z direction
  double * tdm_l = NULL;
  double * tdm_u = NULL;
  tdm.get_l(*tdm_info, &tdm_l);
  tdm.get_u(*tdm_info, &tdm_u);
  const double hz = domain->hz;
  for(size_t k = 0; k < tdm_sizes[0]; k++){
    tdm_l[k] = 1. / hz / hz;
    tdm_u[k] = 1. / hz / hz;
  }
#endif
  return 0;
}

static int init_pencil_rotations(
    const domain_t * domain,
    poisson_solver_t * poisson_solver
){
  const sdecomp_info_t * info = domain->info;
  const size_t r_dsize = sizeof(double);
#if NDIMS == 3
  const size_t c_dsize = sizeof(fftw_complex);
#endif
  if(0 != sdecomp.transpose.construct(info, SDECOMP_X1PENCIL, SDECOMP_Y1PENCIL, r_gl_sizes, r_dsize, &poisson_solver->r_transposer_x1_to_y1)){
    report_failure("SDECOMP x1 to y1 for real");
    return 1;
  }
  if(0 != sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_X1PENCIL, r_gl_sizes, r_dsize, &poisson_solver->r_transposer_y1_to_x1)){
    report_failure("SDECOMP y1 to x1 for real");
    return 1;
  }
#if NDIMS == 3
  if(0 != sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_Z1PENCIL, c_gl_sizes, c_dsize, &poisson_solver->c_transposer_y1_to_z1)){
    report_failure("SDECOMP y1 to z1 for complex");
    return 1;
  }
  if(0 != sdecomp.transpose.construct(info, SDECOMP_Z1PENCIL, SDECOMP_Y1PENCIL, c_gl_sizes, c_dsize, &poisson_solver->c_transposer_z1_to_y1)){
    report_failure("SDECOMP z1 to y1 for complex");
    return 1;
  }
#endif
  return 0;
}

static int init_ffts(
    poisson_solver_t * poisson_solver
){
  const unsigned flags = FFTW_PATIENT | FFTW_DESTROY_INPUT;
  // NOTE: two buffers should be properly given
  //   see "allocate_buffers" above
  // x, real to real
  {
    const int signal_length = r_x1pncl_sizes[SDECOMP_XDIR];
#if NDIMS == 2
    const int repeat_for = r_x1pncl_sizes[SDECOMP_YDIR];
#else
    const int repeat_for = r_x1pncl_sizes[SDECOMP_YDIR] * r_x1pncl_sizes[SDECOMP_ZDIR];
#endif
    fftw_plan * fplan = &poisson_solver->fftw_plan_x[0];
    fftw_plan * bplan = &poisson_solver->fftw_plan_x[1];
    *fplan = fftw_plan_many_r2r(
        1, &signal_length, repeat_for,
        poisson_solver->buf0, NULL, 1, signal_length,
        poisson_solver->buf1, NULL, 1, signal_length,
        (fftw_r2r_kind [1]){FFTW_REDFT10}, flags
    );
    *bplan = fftw_plan_many_r2r(
        1, &signal_length, repeat_for,
        poisson_solver->buf1, NULL, 1, signal_length,
        poisson_solver->buf0, NULL, 1, signal_length,
        (fftw_r2r_kind [1]){FFTW_REDFT01}, flags
    );
    if(NULL == *fplan){
      report_failure("FFTW x-forward");
      return 1;
    }
    if(NULL == *bplan){
      report_failure("FFTW x-backward");
      return 1;
    }
  }
#if NDIMS == 3
  // y, real / complex
  {
    fftw_plan * fplan = &poisson_solver->fftw_plan_y[0];
    fftw_plan * bplan = &poisson_solver->fftw_plan_y[1];
    const int r_signal_length = r_y1pncl_sizes[SDECOMP_YDIR];
    const int c_signal_length = c_y1pncl_sizes[SDECOMP_YDIR];
    const int repeat_for = r_y1pncl_sizes[SDECOMP_ZDIR] * r_y1pncl_sizes[SDECOMP_XDIR];
    *fplan = fftw_plan_many_dft_r2c(
        1, &r_signal_length, repeat_for,
        poisson_solver->buf0, NULL, 1, r_signal_length,
        poisson_solver->buf1, NULL, 1, c_signal_length,
        flags
    );
    *bplan = fftw_plan_many_dft_c2r(
        1, &r_signal_length, repeat_for,
        poisson_solver->buf1, NULL, 1, c_signal_length,
        poisson_solver->buf0, NULL, 1, r_signal_length,
        flags
    );
    if(NULL == *fplan){
      report_failure("FFTW y-forward");
      return 1;
    }
    if(NULL == *bplan){
      report_failure("FFTW y-backward");
      return 1;
    }
  }
#endif
  return 0;
}

static int init_eigenvalues(
    const domain_t * domain,
    poisson_solver_t * poisson_solver
){
  const sdecomp_info_t * info = domain->info;
  double ** evals = &poisson_solver->evals;
#if NDIMS == 2
  // y1 pencil, DCT in x
  const sdecomp_pencil_t pencil = SDECOMP_Y1PENCIL;
  const double signal_lengths[NDIMS - 1] = {
    2. * r_gl_sizes[SDECOMP_XDIR],
  };
  size_t mysizes[NDIMS - 1] = {0};
  sdecomp.get_pencil_mysize(info, pencil, SDECOMP_XDIR, r_gl_sizes[SDECOMP_XDIR], mysizes);
  size_t offsets[NDIMS - 1] = {0};
  sdecomp.get_pencil_offset(info, pencil, SDECOMP_XDIR, r_gl_sizes[SDECOMP_XDIR], offsets);
  const double gridsizes[NDIMS - 1] = {
    domain->lengths[SDECOMP_XDIR] / r_gl_sizes[SDECOMP_XDIR],
  };
  // initialise eigenvalues in homogeneous directions
  *evals = memory_calloc(mysizes[0], sizeof(double));
  for(size_t cnt = 0, i = offsets[0]; i < mysizes[0] + offsets[0]; i++, cnt++){
    (*evals)[cnt] =
      - 4. / pow(gridsizes[0], 2.) * pow(
        sin( g_pi * i / signal_lengths[0] ),
        2.
    );
  }
#else
  // z1 pencil, DCT in x and DFT in y
  const sdecomp_pencil_t pencil = SDECOMP_Z1PENCIL;
  const double signal_lengths[NDIMS - 1] = {
    2. * r_gl_sizes[SDECOMP_XDIR],
    1. * r_gl_sizes[SDECOMP_YDIR],
  };
  size_t mysizes[NDIMS - 1] = {0};
  sdecomp.get_pencil_mysize(info, pencil, SDECOMP_XDIR, c_gl_sizes[SDECOMP_XDIR], mysizes + 0);
  sdecomp.get_pencil_mysize(info, pencil, SDECOMP_YDIR, c_gl_sizes[SDECOMP_YDIR], mysizes + 1);
  size_t offsets[NDIMS - 1] = {0};
  sdecomp.get_pencil_offset(info, pencil, SDECOMP_XDIR, c_gl_sizes[SDECOMP_XDIR], offsets + 0);
  sdecomp.get_pencil_offset(info, pencil, SDECOMP_YDIR, c_gl_sizes[SDECOMP_YDIR], offsets + 1);
  const double gridsizes[NDIMS - 1] = {
    domain->lengths[SDECOMP_XDIR] / r_gl_sizes[SDECOMP_XDIR],
    domain->lengths[SDECOMP_YDIR] / r_gl_sizes[SDECOMP_YDIR],
  };
  // initialise eigenvalues in homogeneous directions
  *evals = memory_calloc(mysizes[0] * mysizes[1], sizeof(double));
  for(size_t cnt = 0, j = offsets[1]; j < mysizes[1] + offsets[1]; j++){
    for(size_t i = offsets[0]; i < mysizes[0] + offsets[0]; i++, cnt++){
      (*evals)[cnt] =
        - 4. / pow(gridsizes[0], 2.) * pow(
          sin( g_pi * i / signal_lengths[0] ),
          2.
        )
        - 4. / pow(gridsizes[1], 2.) * pow(
          sin( g_pi * j / signal_lengths[1] ),
          2.
        );
    }
  }
#endif
  return 0;
}

static int init_poisson_solver(
    const domain_t * domain,
    poisson_solver_t * poisson_solver
){
  // check domain size (global, local, pencils)
  if(0 != compute_pencil_sizes(domain)) return 1;
  // initialise each part of poisson_solver_t
  if(0 != allocate_buffers(poisson_solver))                 return 1;
  if(0 != init_tri_diagonal_solver(domain, poisson_solver)) return 1;
  if(0 != init_pencil_rotations(domain, poisson_solver))    return 1;
  if(0 != init_ffts(poisson_solver))                        return 1;
  if(0 != init_eigenvalues(domain, poisson_solver))         return 1;
  poisson_solver->is_initialised = true;
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(root == myrank){
    printf("DCT-based solver is used\n");
  }
  return 0;
}

static int assign_input(
    const domain_t * domain,
    const double dt_old,
    const double dt_new,
    const fluid_t * fluid,
    double * restrict rhs
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  const double hy = domain->hy;
#if NDIMS == 3
  const double hz = domain->hz;
#endif
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
#if NDIMS == 3
  const double * restrict uz = fluid->uz.data;
#endif
  const double * restrict psi = fluid->psi[0].data;
  const double * restrict den = fluid->den[1].data;
  // normalise FFT beforehand
#if NDIMS == 2
  const double norm = 2. * domain->glsizes[0];
#else
  const double norm = 2. * domain->glsizes[0] * domain->glsizes[1];
#endif
  const double refden = fluid->refden;
  // coefficients in front of each contribution,
  //   0: divergence
  //   1: potential
  // and FFT normalisation
  const double coefs[2] = {
    // rho_ref / gamma dt
    1. / norm / dt_new * refden,
    // c_{n-1}
    1. / norm * -1. / dt_new * dt_old,
  };
#if NDIMS == 2
  for(int cnt = 0, j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++, cnt++){
      // compute local divergence | 13
      const double hx_xm = HXXF(i  );
      const double hx_xp = HXXF(i+1);
      const double jd_xm = JDXF(i  );
      const double jd_x0 = JDXC(i  );
      const double jd_xp = JDXF(i+1);
      const double ux_xm = UX(i  , j  );
      const double ux_xp = UX(i+1, j  );
      const double uy_ym = UY(i  , j  );
      const double uy_yp = UY(i  , j+1);
      const double div = 1. / jd_x0 * (
          - jd_xm / hx_xm * ux_xm + jd_xp / hx_xp * ux_xp
          - jd_x0 / hy    * uy_ym + jd_x0 / hy    * uy_yp
      );
      // additional contribution | 24
      const double den_xm = + 0.5 * DEN(i-1, j  )
                            + 0.5 * DEN(i  , j  );
      const double den_xp = + 0.5 * DEN(i  , j  )
                            + 0.5 * DEN(i+1, j  );
      const double den_ym = + 0.5 * DEN(i  , j-1)
                            + 0.5 * DEN(i  , j  );
      const double den_yp = + 0.5 * DEN(i  , j  )
                            + 0.5 * DEN(i  , j+1);
      const double dpsi_xm = - PSI(i-1, j  )
                             + PSI(i  , j  );
      const double dpsi_xp = - PSI(i  , j  )
                             + PSI(i+1, j  );
      const double dpsi_ym = - PSI(i  , j-1)
                             + PSI(i  , j  );
      const double dpsi_yp = - PSI(i  , j  )
                             + PSI(i  , j+1);
      const double gp_xm = (refden / den_xm - 1.) * dpsi_xm / hx_xm;
      const double gp_xp = (refden / den_xp - 1.) * dpsi_xp / hx_xp;
      const double gp_ym = (refden / den_ym - 1.) * dpsi_ym / hy;
      const double gp_yp = (refden / den_yp - 1.) * dpsi_yp / hy;
      const double add = 1. / jd_x0 * (
          - jd_xm / hx_xm * gp_xm + jd_xp / hx_xp * gp_xp
          - jd_x0 / hy    * gp_ym + jd_x0 / hy    * gp_yp
      );
      rhs[cnt] = (
        + coefs[0] * div
        + coefs[1] * add
      );
    }
  }
#else
  for(int cnt = 0, k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++, cnt++){
        // compute local divergence | 16
        const double hx_xm = HXXF(i  );
        const double hx_xp = HXXF(i+1);
        const double jd_xm = JDXF(i  );
        const double jd_x0 = JDXC(i  );
        const double jd_xp = JDXF(i+1);
        const double ux_xm = UX(i  , j  , k  );
        const double ux_xp = UX(i+1, j  , k  );
        const double uy_ym = UY(i  , j  , k  );
        const double uy_yp = UY(i  , j+1, k  );
        const double uz_zm = UZ(i  , j  , k  );
        const double uz_zp = UZ(i  , j  , k+1);
        const double div = 1. / jd_x0 * (
            - jd_xm / hx_xm * ux_xm + jd_xp / hx_xp * ux_xp
            - jd_x0 / hy    * uy_ym + jd_x0 / hy    * uy_yp
            - jd_x0 / hz    * uz_zm + jd_x0 / hz    * uz_zp
        );
        // additional contribution | 35
        const double den_xm = + 0.5 * DEN(i-1, j  , k  )
                              + 0.5 * DEN(i  , j  , k  );
        const double den_xp = + 0.5 * DEN(i  , j  , k  )
                              + 0.5 * DEN(i+1, j  , k  );
        const double den_ym = + 0.5 * DEN(i  , j-1, k  )
                              + 0.5 * DEN(i  , j  , k  );
        const double den_yp = + 0.5 * DEN(i  , j  , k  )
                              + 0.5 * DEN(i  , j+1, k  );
        const double den_zm = + 0.5 * DEN(i  , j  , k-1)
                              + 0.5 * DEN(i  , j  , k  );
        const double den_zp = + 0.5 * DEN(i  , j  , k  )
                              + 0.5 * DEN(i  , j  , k+1);
        const double dpsi_xm = - PSI(i-1, j  , k  )
                               + PSI(i  , j  , k  );
        const double dpsi_xp = - PSI(i  , j  , k  )
                               + PSI(i+1, j  , k  );
        const double dpsi_ym = - PSI(i  , j-1, k  )
                               + PSI(i  , j  , k  );
        const double dpsi_yp = - PSI(i  , j  , k  )
                               + PSI(i  , j+1, k  );
        const double dpsi_zm = - PSI(i  , j  , k-1)
                               + PSI(i  , j  , k  );
        const double dpsi_zp = - PSI(i  , j  , k  )
                               + PSI(i  , j  , k+1);
        const double gp_xm = (refden / den_xm - 1.) * dpsi_xm / hx_xm;
        const double gp_xp = (refden / den_xp - 1.) * dpsi_xp / hx_xp;
        const double gp_ym = (refden / den_ym - 1.) * dpsi_ym / hy;
        const double gp_yp = (refden / den_yp - 1.) * dpsi_yp / hy;
        const double gp_zm = (refden / den_zm - 1.) * dpsi_zm / hz;
        const double gp_zp = (refden / den_zp - 1.) * dpsi_zp / hz;
        const double add = 1. / jd_x0 * (
            - jd_xm / hx_xm * gp_xm + jd_xp / hx_xp * gp_xp
            - jd_x0 / hy    * gp_ym + jd_x0 / hy    * gp_yp
            - jd_x0 / hz    * gp_zm + jd_x0 / hz    * gp_zp
        );
        rhs[cnt] = (
          + coefs[0] * div
          + coefs[1] * add
        );
      }
    }
  }
#endif
  return 0;
}

static int extract_output(
    const domain_t * domain,
    const double * restrict rhs,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  double * restrict psi = fluid->psi[1].data;
#if NDIMS == 2
  for(int cnt = 0, j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++, cnt++){
      PSI(i, j) = rhs[cnt];
    }
  }
#else
  for(int cnt = 0, k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++, cnt++){
        PSI(i, j, k) = rhs[cnt];
      }
    }
  }
#endif
  fluid_update_boundaries_psi(domain, &fluid->psi[1]);
  return 0;
}

static int solve_linear_systems(
    poisson_solver_t * poisson_solver
){
  // size of system (length) and how many such systems to be solved
  // NOTE: although size_of_system is the same as tdm.get_size gives,
  //   repeat_for is different from what tdm.get_nrhs returns (=1)
  //   here repeat_for is the degree of freedom in the wavespace
  const size_t size_of_system = poisson_solver->tdm_sizes[0];
  const size_t repeat_for     = poisson_solver->tdm_sizes[1];
  // tri-diagonal matrix
  tdm_info_t * tdm_info = poisson_solver->tdm_info;
  double * restrict tdm_l = NULL;
  double * restrict tdm_u = NULL;
  double * restrict tdm_c = NULL;
  tdm.get_l(tdm_info, &tdm_l);
  tdm.get_u(tdm_info, &tdm_u);
  tdm.get_c(tdm_info, &tdm_c);
  // eigenvalues coming from Fourier projection
  const double * restrict evals = poisson_solver->evals;
#if NDIMS == 2
  double * restrict rhs = poisson_solver->buf0;
#else
  fftw_complex * restrict rhs = poisson_solver->buf0;
#endif
  for(size_t m = 0; m < repeat_for; m++){
    // set center diagonal components
    for(size_t n = 0; n < size_of_system; n++){
      tdm_c[n] = - tdm_l[n] - tdm_u[n] + evals[m];
    }
    tdm.solve(tdm_info, rhs + m * size_of_system);
  }
  return 0;
}

/**
 * @brief compute scalar potential psi to correct velocity
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     dt_old : previous time step size
 * @param[in]     dt_new : current  time step size
 * @param[in,out] fluid  : velocity (in), scalar potential psi (out)
 * @return               : (success) 0
 *                       : (failure) 1
 */
int fluid_compute_potential(
    const domain_t * domain,
    const double dt_old,
    const double dt_new,
    fluid_t * fluid
){
  static poisson_solver_t poisson_solver = {
    .is_initialised = false,
  };
  // initialise Poisson solver
  if(!poisson_solver.is_initialised){
    if(0 != init_poisson_solver(domain, &poisson_solver)){
      // failed to initialise Poisson solver
      return 1;
    }
  }
  // compute right-hand side of Poisson equation
  // assigned to buf0
  assign_input(domain, dt_old, dt_new, fluid, poisson_solver.buf0);
  // solve the equation
  // project x to wave space
  // f(x, y)    -> f(k_x, y)
  // f(x, y, z) -> f(k_x, y, z)
  // from buf0 to buf1
  fftw_execute(poisson_solver.fftw_plan_x[0]);
  // transpose real x1pencil to y1pencil
  // from buf1 to buf0
  sdecomp.transpose.execute(
      poisson_solver.r_transposer_x1_to_y1,
      poisson_solver.buf1,
      poisson_solver.buf0
  );
#if NDIMS == 3
  // project y to wave space
  // f(k_x, y, z) -> f(k_x, k_y, z)
  // from buf0 to buf1
  fftw_execute(poisson_solver.fftw_plan_y[0]);
  // transpose complex y1pencil to z1pencil
  // from buf1 to buf0
  sdecomp.transpose.execute(
      poisson_solver.c_transposer_y1_to_z1,
      poisson_solver.buf1,
      poisson_solver.buf0
  );
#endif
  // solve linear systems
  solve_linear_systems(&poisson_solver);
#if NDIMS == 3
  // transpose complex z1pencil to y1pencil
  // from buf0 to buf1
  sdecomp.transpose.execute(
      poisson_solver.c_transposer_z1_to_y1,
      poisson_solver.buf0,
      poisson_solver.buf1
  );
  // project y to physical space
  // f(k_x, k_y, z) -> f(k_x, y, z)
  // from buf1 to buf0
  fftw_execute(poisson_solver.fftw_plan_y[1]);
#endif
  // transpose real y1pencil to x1pencil
  // from buf0 to buf1
  sdecomp.transpose.execute(
      poisson_solver.r_transposer_y1_to_x1,
      poisson_solver.buf0,
      poisson_solver.buf1
  );
  // project x to physical space
  // f(k_x, y)    -> f(x, y)
  // f(k_x, y, z) -> f(x, y, z)
  // from buf1 to buf0
  fftw_execute(poisson_solver.fftw_plan_x[1]);
  extract_output(domain, poisson_solver.buf0, fluid);
  return 0;
}

