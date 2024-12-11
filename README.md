# Numerical Methods in Astrophysics


Here you will find Python projects which have been inspired by the book "Numerical Methods in Astrophysics" by P. Bodenheimer, et al.

Each project will cover a different numerical method.


## Internal Viscosity (Numerical Diffusion)

This simple solution (1st order) to the 1D continuity equation demonstrates the effect of internal viscosity spreading an initial discontinuity into smooth, broad ramps as time passes. (In the figure, v = 100 and ∆t = 0.002 = one fifth the stability limit of ∆x/v)!


## N-Body Problem (Euler)

This script solves the N-body problem in 3D using Euler's method by repeatedly applying a finite difference approximation. This serves as a simple incarnation of an N-body algorithm. The values used in the script, and thus to produce the figure, are for the solar system evolving over 1 Earth-year. I'm only considering the Sun, the Earth, and Mars; thus, a 3-body problem. However, this program can be simply altered to include more bodies at different starting positions and trajectories, making it applicable for more complicated N-body (small N) problems!


## 1D Lagrangian Hydrodynamics (Shock)

This project tackles the explosion of a supernova into the ISM. It solves, on a Lagrangian grid, the relevant 1D adiabatic hydrodynamic equations (including artifical viscosity). There are 4 test cases (given as inputs in lh1.txt) which produce 4 of the figures; different shocks are shown by these varying input parameters, all solved using the lh1.py script. The fifth figure, "nx = 500.png", demonstrates the features becoming very diffused with poor energy conservation due to the small number of grid points. (Finally, lh1.out and lh2.out are just output data of the script which can be read for plotting).


## More to come...
