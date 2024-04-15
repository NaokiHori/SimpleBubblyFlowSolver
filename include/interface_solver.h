#if !defined(INTERFACE_SOLVER_H)
#define INTERFACE_SOLVER_H

#include "domain.h"
#include "fluid.h"
#include "interface.h"

extern int interface_init(
    const char dirname_ic[],
    const domain_t * domain,
    interface_t * interface
);

extern int interface_save(
    const char dirname[],
    const domain_t * domain,
    const interface_t * interface
);

extern int interface_compute_curvature_tensor(
    const domain_t * domain,
    interface_t * interface
);

extern int interface_compute_force(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
);

extern int interface_compute_mass_flux(
    const domain_t * domain,
    const fluid_t * fluid,
    interface_t * interface
);

extern int interface_update_vof(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    const fluid_t * fluid,
    interface_t * interface
);

extern int interface_update_boundaries_vof(
    const domain_t * domain,
    array_t * vof
);

#endif // INTERFACE_SOLVER_H
