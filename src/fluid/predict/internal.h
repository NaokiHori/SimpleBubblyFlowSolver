#if !defined(FLUID_COMPUTE_RHS_INTERNAL)
#define FLUID_COMPUTE_RHS_INTERNAL

#include "interface.h"

extern int compute_txx(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_txy(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_txz(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_tyy(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_tyz(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_tzz(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_rhs_ux(
    const domain_t * domain,
    fluid_t * fluid,
    const interface_t * interface
);

extern int compute_rhs_uy(
    const domain_t * domain,
    fluid_t * fluid,
    const interface_t * interface
);

extern int compute_rhs_uz(
    const domain_t * domain,
    fluid_t * fluid,
    const interface_t * interface
);

extern int predict_ux(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

extern int predict_uy(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

extern int predict_uz(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

#endif // FLUID_COMPUTE_RHS_INTERNAL
