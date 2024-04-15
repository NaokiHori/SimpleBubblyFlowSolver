#if !defined(RUNGE_KUTTA_H)
#define RUNGE_KUTTA_H

#include <stdint.h>

// Runge-Kutta configurations
// indices
extern const uint_fast8_t rk_a; // 0
extern const uint_fast8_t rk_b; // 1
extern const uint_fast8_t rk_g; // 2
// coefficients of three-step RK scheme
// NOTE: only three is allowed
#define RKSTEPMAX 3
typedef double rkcoef_t[RKSTEPMAX];
extern const rkcoef_t rkcoefs[RKSTEPMAX];

#endif // RUNGE_KUTTA_H
