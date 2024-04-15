#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "interface.h"
#include "interface_solver.h"
#include "integrate.h"
#include "decide_dt.h"

// integrate the equations for one time step
int integrate(
    const domain_t * domain,
    fluid_t * fluid,
    interface_t * interface,
    double * dt
){
  // decide time step size
  if(0 != decide_dt(domain, fluid, interface, dt)){
    return 1;
  }
  // Runge-Kutta iterations
  // max iteration, should be three
  for(size_t rkstep = 0; rkstep < RKSTEPMAX; rkstep++){
    // compute k-step density and viscosity
    if(0 != fluid_compute_density(domain, fluid, interface, 0)){
      return 1;
    }
    if(0 != fluid_compute_viscosity(domain, fluid, interface)){
      return 1;
    }
    // update vof field
    if(0 != interface_compute_curvature_tensor(domain, interface)){
      return 1;
    }
    if(0 != interface_compute_force(domain, fluid, interface)){
      return 1;
    }
    if(0 != interface_update_vof(domain, rkstep, *dt, fluid, interface)){
      return 1;
    }
    // compute (k+1)-step density
    if(0 != fluid_compute_density(domain, fluid, interface, 1)){
      return 1;
    }
    // override vof flux by mass flux
    if(0 != interface_compute_mass_flux(domain, fluid, interface)){
      return 1;
    }
    // predict flow field
    if(0 != fluid_predict_field(domain, rkstep, *dt, fluid, interface)){
      return 1;
    }
    // local time step sizes (gamma dt)
    const double dt_old = rkcoefs[(rkstep + 2) % 3][rk_g] * (*dt);
    const double dt_new = rkcoefs[(rkstep + 0) % 3][rk_g] * (*dt);
    // compute scalar potential
    if(0 != fluid_compute_potential(domain, dt_old, dt_new, fluid)){
      return 1;
    }
    // correct velocity field to satisfy mass conservation
    if(0 != fluid_correct_velocity(domain, dt_old, dt_new, fluid)){
      return 1;
    }
    // update pressure
    if(0 != fluid_update_pressure(domain, fluid)){
      return 1;
    }
  }
  return 0;
}

