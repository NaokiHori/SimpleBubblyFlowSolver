#include "sdecomp.h"
#include "array.h"
#include "domain.h"
#include "interface.h"
#include "fileio.h"

int interface_save(
    const char dirname[],
    const domain_t * domain,
    const interface_t * interface
){
  // serial
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(root == myrank){
    fileio.w_serial(dirname, "We", 0, NULL, fileio.npy_double, sizeof(double), &interface->We);
  }
  // collective
  array.dump(domain, dirname, "vof", fileio.npy_double, &interface->vof);
  return 0;
}

