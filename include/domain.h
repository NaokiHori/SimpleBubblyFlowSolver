#if !defined(DOMAIN_H)
#define DOMAIN_H

#include "sdecomp.h"

// definition of a structure domain_t
/**
 * @struct domain_t
 * @brief struct storing parameters relevant to spatial domain
 * @var info       : MPI domain decomposition
 * @var glsizes    : global     number of grid points in each direction
 * @var mysizes    : local (my) number of grid points in each direction
 * @var offsets    : offsets to my starting index in each direction
 * @var lengths    : domain size in each direction
 * @var xf, xc     : cell-face and cell-center locations in x direction
 * @var jdxf, jdxc : Jacobian determinants at x cell faces and centers
 * @var hxxf, hxxc : scale factors at x cell faces and centers
 * @var hy, hz     : y and z scale factors
 */
typedef struct {
  sdecomp_info_t * info;
  size_t glsizes[NDIMS];
  size_t mysizes[NDIMS];
  size_t offsets[NDIMS];
  double lengths[NDIMS];
  double * restrict xf, * restrict xc;
  double * restrict jdxf, * restrict jdxc;
  double * restrict hxxf, * restrict hxxc;
  double hy;
} domain_t;

// constructor
extern int domain_init(
    const char dirname_ic[],
    domain_t * domain
);

// save members which are necessary to restart
extern int domain_save(
    const char dirname[],
    const domain_t * domain
);

#endif // DOMAIN_H
