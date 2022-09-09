# TaylorCouette-AFiD
Clean version of Taylor-Couette code with AFiD routine names

There are three velocity components, azimuthal-radial-axial. The radial velocity component Q_r is written as r*v_r as proposed in Verzicco and Orlandi (JCP, 1996). The code is written in a semiimplicit manner, and paralelized according to the same directives as the Rayleigh-Benard AFiD code (van der Poel et al., C&F 2015).

The code runs by reading the options in bou.in: 

NTM/NZM/NRM, resolution in the theta-z-r direction

NSST - Choose between 3rd order Runge-Kutta (3) or second order Adams-Bashworth (1)

IDTV - Variable timestep (1) or fixed timestep (0)

NREAD - 0 means restart simulation, 1 means read old files

NTST - maximum number of timesteps before quitting

TOUT - simulation time at which the code outputs data & calculates statistics

TMAX - maximum simulation time before quitting

IRESET - reset logs if 1, if 0 append

WALLTIMEMAX - maximum wall time (real world time) before quitting

The code quits on either NTST, TMAX or WALLTIMEMAX, whichever of the three is satisfied before

ALX3D - Periodic length in the axial direction

RINT - Inner cylinder radius, leave as 1

REXT - Outer cylinder radius

ISTR/STR3 - Grid clustering near the inner/outer cylinders

LAMB - Rotational symmetry imposed

VIROTINT/VIROTEXT - rotation velocity of inner and outer cylinder

AMPL_PERT - amplitude of the inner cylinder oscillation

FREQ_PERT - frequency of the inner cylinder perturbation (in time, not in rad/s)

The inner cylinder velocity is VIROTINT + AMPL_PERT · sin(2pi · FREQ_PERT · time)

1/NU - inverse of the kinematic viscosity

ROTP - rotation parameter

RESID - maximum divergence of the field before crashing

CFLMAX - maximum cfl to determine time-step in variable time-step

TSTA - time at which the code begins to take statistics

STAREAD - read the statistics from the previous run, and append (the code outputs in stats the mean & rms velocity fields)

DTMIN - minimum/maximum time step size for variable time-step

VMAX - maximum velocity before code stop

CFLLIM - maximum CFL for fixed time-step before code stop

Output options are the following:

1. Two-dimensional slabs, outputted every TFRAME2D

(a) Constant R cuts are enabled by the NOUT_RCUT flag in bou.in, which sets the number of cuts to output, from 0 to 9. The approximate location of the rcuts are set in rcuts.in, and the real r-index of the cuts is outputted in rcuts/rcutslocs.out. 

(b) Constant Azimuth cuts are enabled by the OUT_THCUT and are outputted in thcuts/

(c) Streamwise (azimuthal) averages are enabled by OUT_STWAVE and are outputted in stave/

2. Full three dimensional flowfields, outputted every TFRAME3D, and are enabled by OUT_3DFIELD.
