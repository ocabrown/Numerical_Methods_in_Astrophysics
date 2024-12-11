# Numerical_Methods_in_Astrophysics


Here you will find Python projects which have been inspired by the book "Numerical Methods in Astrophysics" by P. Bodenheimer, et al.

Each project will cover a different numerical method.


## Internal Viscosity (Numerical Diffusion)

This simple solution (1st order) to the 1D continuity equation demonstrates the effect of internal viscosity spreading an initial discontinuity into smooth, broad ramps as time passes. (In the figure, v = 100 and ∆t = 0.002 = one fifth the stability limit of ∆x/v)!


## N-Body Problem (Euler)

This script solves the N-body problem in 3D using Euler's method by repeatedly applying a finite difference approximation. This serves as a simple incarnation of an N-body algorithm. The values used in the script, and thus to produce the figure, are for the solar system evolving over 1 Earth-year. I'm only considering the Sun, the Earth, and Mars; thus, a 3-body problem. However, this program can be simply altered to include more bodies at different starting positions and trajectories, making it applicable for more complicated N-body (small N) problems!


