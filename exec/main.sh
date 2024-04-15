#!/bin/bash

## temporal information
# maximum duration (in free-fall time)
export timemax=5.0e+0
# maximum duration (in wall time [s])
export wtimemax=6.0e+2
# logging rate (in free-fall time)
export log_rate=1.0e-1
# save rate (in free-fall time)
export save_rate=5.0e-1
# save after (in free-fall time)
export save_after=0.0e+0
# statistics collection rate (in free-fall time)
export stat_rate=1.0e-1
# statistics collection after (in free-fall time)
export stat_after=1.0e+2

## safety factors to decide time step size
## for advective and diffusive terms
export coef_dt_adv=0.20
export coef_dt_dif=0.50
export coef_dt_int=0.95

## physical parameters
export Re=35.
export We=10.
export Fr=1.
export denr=0.0013
export visr=0.018

# give name of the directory in which the initial conditions
#   (incl. domain size etc.) are stored as an argument
dirname_ic=initial_condition/output
# dirname_ic=$(find output/save -type d | sort | tail -n 1)

mpirun -n 2 --oversubscribe ./a.out ${dirname_ic}
