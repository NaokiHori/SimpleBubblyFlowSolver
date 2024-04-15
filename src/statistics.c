#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "sdecomp.h"
#include "param.h"
#include "memory.h"
#include "domain.h"
#include "fluid.h"
#include "interface.h"
#include "statistics.h"
#include "fileio.h"
#include "config.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#if NDIMS == 3
#include "array_macros/fluid/uz.h"
#endif
#include "array_macros/interface/vof.h"
#include "array_macros/statistics/ux1.h"
#include "array_macros/statistics/uy1.h"
#if NDIMS == 3
#include "array_macros/statistics/uz1.h"
#endif
#include "array_macros/statistics/vof1.h"

// parameters to specify directory name
static const char g_dirname_prefix[] = {"output/stat/step"};
static const int g_dirname_ndigits = 10;

// name of directory
static char * g_dirname = NULL;
static size_t g_dirname_nchars = 0;

// scheduler
static double g_rate = 0.;
static double g_next = 0.;

// data
static size_t g_num = 0;
static array_t g_ux1 = {0};
static array_t g_uy1 = {0};
#if NDIMS == 3
static array_t g_uz1 = {0};
#endif
static array_t g_vof1 = {0};

/**
 * @brief constructor - initialise and allocate internal buffers, schedule collection
 * @param[in] domain : information about domain decomposition and size
 * @param[in] time   : current time (hereafter in free-fall time units)
 * @return           : error code
 */
static int init(
    const domain_t * domain,
    const double time
){
  // fetch timings
  if(0 != config.get_double("stat_rate", &g_rate)){
    return 1;
  }
  double after = 0.;
  if(0 != config.get_double("stat_after", &after)){
    return 1;
  }
  g_next = g_rate * ceil(
      fmax(DBL_EPSILON, fmax(time, after)) / g_rate
  );
  // allocate directory name
  g_dirname_nchars =
    + strlen(g_dirname_prefix)
    + g_dirname_ndigits;
  g_dirname = memory_calloc(g_dirname_nchars + 2, sizeof(char));
  // prepare arrays
  if(0 != array.prepare(domain, UX1_NADDS,  sizeof(double), &g_ux1))  return 1;
  if(0 != array.prepare(domain, UY1_NADDS,  sizeof(double), &g_uy1))  return 1;
#if NDIMS == 3
  if(0 != array.prepare(domain, UZ1_NADDS,  sizeof(double), &g_uz1))  return 1;
#endif
  if(0 != array.prepare(domain, VOF1_NADDS, sizeof(double), &g_vof1)) return 1;
  // report
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(root == myrank){
    FILE * stream = stdout;
    fprintf(stream, "STATISTICS\n");
    fprintf(stream, "\tdest: %s\n", g_dirname_prefix);
    fprintf(stream, "\tnext: % .3e\n", g_next);
    fprintf(stream, "\trate: % .3e\n", g_rate);
    fflush(stream);
  }
  return 0;
}

/**
 * @brief getter of a member: g_next
 * @return : g_next
 */
static double get_next_time(
    void
){
  return g_next;
}

/**
 * @brief compute ux^1 and add results to the array
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] ux     : x velocity
 */
static void collect_mean_ux(
    const domain_t * domain,
    const double * restrict ux
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  double * restrict ux1 = g_ux1.data;
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize + 1; i++){
      UX1(i, j) += UX(i, j);
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize + 1; i++){
        UX1(i, j, k) += UX(i, j, k);
      }
    }
  }
#endif
}

/**
 * @brief compute uy^1 and add results to the array
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] uy     : y velocity
 */
static void collect_mean_uy(
    const domain_t * domain,
    const double * restrict uy
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  double * restrict uy1 = g_uy1.data;
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 0; i <= isize + 1; i++){
      UY1(i, j) += UY(i, j);
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize + 1; i++){
        UY1(i, j, k) += UY(i, j, k);
      }
    }
  }
#endif
}

#if NDIMS == 3
/**
 * @brief compute uz^1 and add results to the array
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] uz     : z velocity
 */
static void collect_mean_uz(
    const domain_t * domain,
    const double * restrict uz
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  double * restrict uz1 = g_uz1.data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize + 1; i++){
        UZ1(i, j, k) += UZ(i, j, k);
      }
    }
  }
}
#endif

static void collect_mean_vof(
    const domain_t * domain,
    const double * restrict vof
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  double * restrict vof1 = g_vof1.data;
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 0; i <= isize + 1; i++){
      VOF1(i, j) += VOF(i, j);
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize + 1; i++){
        VOF1(i, j, k) += VOF(i, j, k);
      }
    }
  }
#endif
}

/**
 * @brief accumulate statistical data
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] fluid  : flow field
 * @return           : error code
 */
static int collect(
    const domain_t * domain,
    const fluid_t * fluid,
    const interface_t * interface
){
  // collect temporally-averaged quantities
  collect_mean_ux(domain, fluid->ux.data);
  collect_mean_uy(domain, fluid->uy.data);
#if NDIMS == 3
  collect_mean_uz(domain, fluid->uz.data);
#endif
  collect_mean_vof(domain, interface->vof.data);
  // increment number of samples
  g_num += 1;
  // schedule next event
  g_next += g_rate;
  return 0;
}

/**
 * @brief save structures which contains collected statistical data
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] step   : current time step
 * @return           : error code
 */
static int output(
    const domain_t * domain,
    const size_t step
){
  // when no statistics are collected (g_num is 0),
  //   no reason to save, so abort
  if(0 == g_num){
    return 0;
  }
  // set directory name
  snprintf(
      g_dirname,
      g_dirname_nchars + 1,
      "%s%0*zu",
      g_dirname_prefix,
      g_dirname_ndigits,
      step
  );
  // get communicator to identify the main process
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  // create directory and save scalars from main process
  if(root == myrank){
    // although it may fail, anyway continue, which is designed to be safe
    fileio.mkdir(g_dirname);
    // save scalars
    fileio.w_serial(g_dirname, "num", 0, NULL, fileio.npy_size_t, sizeof(size_t), &g_num);
  }
  // wait for the main process to complete making directory
  MPI_Barrier(MPI_COMM_WORLD);
  // save domain info (coordinates)
  domain_save(g_dirname, domain);
  // save collected statistics
  typedef struct {
    const char * name;
    const array_t * array;
  } variable_t;
  const variable_t variables[] = {
    {.name = "ux1",  .array = &g_ux1},
    {.name = "uy1",  .array = &g_uy1},
#if NDIMS == 3
    {.name = "uz1",  .array = &g_uz1},
#endif
    {.name = "vof1", .array = &g_vof1},
  };
  for(size_t index = 0; index < sizeof(variables) / sizeof(variable_t); index++){
    const variable_t * v = variables + index;
    array.dump(domain, g_dirname, v->name, fileio.npy_double, v->array);
  }
  return 0;
}

const statistics_t statistics = {
  .init          = init,
  .collect       = collect,
  .output        = output,
  .get_next_time = get_next_time,
};

