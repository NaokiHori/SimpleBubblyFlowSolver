#if !defined(FLUID_H)
#define FLUID_H

#include "array.h"
#include "domain.h"

// definition of a structure fluid_t_
/**
 * @struct fluid_t
 * @brief struct storing fluid-related variables
 * @var ux, uy, uz : velocity in each direction
 * @var p, psi     : pressure, scalar potentials
 * @var den        : density
 * @var visux      : dynamic viscosity at ux
 * @var visuy      : dynamic viscosity at uy
 * @var visuz      : dynamic viscosity at uz
 * @var srcux      : Runge-Kutta source terms for ux
 * @var srcuy      : Runge-Kutta source terms for uy
 * @var srcuz      : Runge-Kutta source terms for uz
 * @var t[x-z][x-z]: shear-stress tensor
 * @var Re         : Reynolds number
 * @var Fr         : Froude   number
 * @var denr       : density ratio
 * @var visr       : viscosity ratio
 * @var refden     : reference density
 */
typedef struct {
  array_t ux;
  array_t uy;
  array_t uz;
  array_t p;
  array_t psi[2];
  array_t den[2];
  array_t visux;
  array_t visuy;
  array_t visuz;
  array_t txx;
  array_t txy;
  array_t txz;
  array_t tyy;
  array_t tyz;
  array_t tzz;
  array_t srcux[3];
  array_t srcuy[3];
  array_t srcuz[3];
  double Re;
  double Fr;
  double denr;
  double visr;
  double refden;
} fluid_t;

#endif // FLUID_H
