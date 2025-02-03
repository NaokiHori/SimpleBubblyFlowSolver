#if !defined(PARAM_H)
#define PARAM_H

// fixed parameters, which are usually fixed
//   but still user can easily control, are declared
// they are defined under src/param/xxx.c

#include <stdbool.h>

/* boundary-condition.c */
// NOTE: impermeable walls and Neumann BC for the pressure are unchangeable
// negative-x-wall velocity in y direction
extern const double param_uy_xm;
// positive-x-wall velocity in y direction
extern const double param_uy_xp;
// negative-x-wall velocity in z direction
extern const double param_uz_xm;
// positive-x-wall velocity in z direction
extern const double param_uz_xp;

#endif // PARAM_H
