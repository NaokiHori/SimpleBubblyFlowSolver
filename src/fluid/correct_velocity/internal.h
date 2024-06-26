#if !defined(FLUID_CORRECT_VELOCITY_INTERNAL_H)
#define FLUID_CORRECT_VELOCITY_INTERNAL_H

extern int fluid_correct_velocity_ux(
    const domain_t * domain,
    const double dt_old,
    const double dt_new,
    fluid_t * fluid
);

extern int fluid_correct_velocity_uy(
    const domain_t * domain,
    const double dt_old,
    const double dt_new,
    fluid_t * fluid
);

#endif // FLUID_CORRECT_VELOCITY_INTERNAL_H
