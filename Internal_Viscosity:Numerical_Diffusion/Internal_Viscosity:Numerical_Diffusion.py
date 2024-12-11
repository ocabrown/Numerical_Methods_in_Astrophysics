import numpy as np
import matplotlib.pyplot as plt

# Parameters
v = 100.0           # Speed
dt = 0.002          # Time step
dx = 1.0            # Space step (modifiable)
x_max = 400         # Maximum x value
t_max = 3.0         # Maximum time value
d = v * dt / dx     # Courant number

# Initialize the grid
x = np.arange(0, x_max + dx, dx)
n_x = len(x)
n_t = int(t_max / dt) + 1   # Total number of time steps

# Create a 2D array to store density at each time step
p = np.zeros((n_t, n_x))

# Initial condition for p(x, 0)
p[0, 0] = 1.0                   # p(0,0) = 1
p[0, 1:int(50/dx)+1] = 0.5      # p(x,0) = 0.5 for 0 < x <= 50
p[0, int(50/dx)+1:] = 0.0       # p(x,0) = 0 for x > 50

# Time evolution
for n in range(1, n_t):
    for j in range(1, n_x - 1):
        # Apply the explicit scheme
        p[n, j] = -d * (p[n-1, j] - p[n-1, j-1]) + p[n-1, j]

    # Update the boundary conditions if needed
    p[n, 0] = 1.0   # p(0, t) remains 1.0
    p[n, -1] = 0.0  # Boundary condition at x_max

# Times to plot
times_to_plot = [0, 1, 2, 3]                        # Times at which we will plot
plot_steps = [int(t / dt) for t in times_to_plot]   # Convert times to steps

# Plotting the results at different time steps
plt.figure(figsize=(10, 6))
for step in plot_steps:
    plt.plot(x, p[step], label=f"t = {step * dt:.1f} s")

plt.title("1D Continuity Equation Evolution")
plt.xlabel("x")
plt.ylabel("Density (p)")
plt.legend()
plt.grid(True)
plt.show()
