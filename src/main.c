#include <stdio.h>
#include <mpi.h>
#include "timer.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "interface.h"
#include "interface_solver.h"
#include "integrate.h"
#include "statistics.h"
#include "save.h"
#include "logging.h"
#include "config.h"
#include "fileio.h"

static int save_entrypoint(
    const domain_t * domain,
    const size_t step,
    const double time,
    const fluid_t * fluid,
    const interface_t * interface
){
  char * dirname = NULL;
  save.prepare(domain, step, &dirname);
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(root == myrank){
    fileio.w_serial(dirname, "step", 0, NULL, fileio.npy_size_t, sizeof(size_t), &step);
    fileio.w_serial(dirname, "time", 0, NULL, fileio.npy_double, sizeof(double), &time);
  }
  domain_save(dirname, domain);
  fluid_save(dirname, domain, fluid);
  interface_save(dirname, domain, interface);
  return 0;
}

/**
 * @brief main function
 * @param[in] argc : number of arguments (expect 2)
 * @param[in] argv : name of the directory
 *                     where a set of initial condition is contained
 * @return         : error code
 */
int main(
    int argc,
    char * argv[]
){
  // launch MPI, start timer
  MPI_Init(NULL, NULL);
  const int root = 0;
  int myrank = root;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  const double tic = timer();
  // find name of directory where IC is stored
  if(2 != argc){
    if(root == myrank){
      printf("directory name should be given as input\n");
    }
    goto abort;
  }
  const char * dirname_ic = argv[1];
  // initialise fileio object
  if(0 != fileio.init()){
    goto abort;
  }
  // initialise time step and time units
  size_t step = 0;
  if(0 != fileio.r_serial(dirname_ic, "step", 0, NULL, fileio.npy_size_t, sizeof(size_t), &step)){
    goto abort;
  }
  double time = 0.;
  if(0 != fileio.r_serial(dirname_ic, "time", 0, NULL, fileio.npy_double, sizeof(double), &time)){
    goto abort;
  }
  // initialise structures
  domain_t domain = {0};
  if(0 != domain_init(dirname_ic, &domain)){
    goto abort;
  }
  fluid_t fluid = {0};
  if(0 != fluid_init(dirname_ic, &domain, &fluid)){
    goto abort;
  }
  interface_t interface = {0};
  if(0 != interface_init(dirname_ic, &domain, &interface)){
    goto abort;
  }
  // initialise auxiliary objects
  if(0 != logging.init(&domain, time)){
    goto abort;
  }
  if(0 != save.init(&domain, time)){
    goto abort;
  }
  if(0 != statistics.init(&domain, time)){
    goto abort;
  }
  // check termination conditions
  double timemax = 0.;
  if(0 != config.get_double("timemax", &timemax)){
    goto abort;
  }
  double wtimemax = 0.;
  if(0 != config.get_double("wtimemax", &wtimemax)){
    goto abort;
  }
  // report
  if(root == myrank){
    printf("step: %zu, time: % .7e\n", step, time);
    printf("timemax: % .7e, wtimemax: % .7e\n", timemax, wtimemax);
  }
  // main loop
  for(;;){
    // time step size
    double dt = 0.;
    // integrate for one time step
    if(0 != integrate(&domain, &fluid, &interface, &dt)){
      goto abort;
    }
    // update step and simulation / wall time
    step += 1;
    time += dt;
    const double toc = timer();
    // terminate if one of the following conditions is met
    // the simulation is finished
    if(timemax < time){
      break;
    }
    // wall time limit is reached
    if(wtimemax < toc - tic){
      break;
    }
    // compute and output log regularly
    if(logging.get_next_time() < time){
      logging.check_and_output(&domain, step, time, dt, toc - tic, &fluid, &interface);
    }
    // save flow fields regularly
    if(save.get_next_time() < time){
      save_entrypoint(&domain, step, time, &fluid, &interface);
    }
    // collect statistics regularly
    if(statistics.get_next_time() < time){
      statistics.collect(&domain, &fluid, &interface);
    }
  }
  // finalisation
  // save final flow fields
  save_entrypoint(&domain, step, time, &fluid, &interface);
  // save collected statistics
  statistics.output(&domain, step);
  // finalise MPI
abort:
  MPI_Finalize();
  return 0;
}

