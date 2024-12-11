import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-8  # Gravitational constant in cgs

# Number of bodies
N = 3

# Initialize masses for the 2 point masses
mass = np.array([1.989e33, 5.972e27, 6.417e26])    # Masses in g (Sun, Earth, Mars)

# Initialize positions (x) and velocities (v) in 3D space
x_vals = np.array([[0, 0, 0], [1.496e13, 0, 0], [2.251e13, 0, 0]])   # Initial Positions in cm
v_vals = np.array([[0, 0, 0], [0, 2.978e6, 0], [0, 2.408e6, 0]])     # Initial Orbital Velocities in cm/s

# Initial state vector W = [x_1, v_1, ..., x_N, v_N]
W_vals = np.hstack((x_vals, v_vals))   # Stack positions and velocities horizontally

# Time step
h_val = 1000    # Time step in seconds


# W = [x, v] is the 6N dimensional vector where W[i] = [x[i], v[i]]
def find_acceleration(m, x):
    """ Compute the gravitational acceleration for each body based on the positions of all bodies. """

    a = np.zeros((N, 3))    # Initialize accelerations for all bodies (3D)

    for i in range(N):
        a_i = np.zeros(3)   # The acceleration acting on body i

        for j in range(N):
            if i != j:
                # Vector from body j to body i
                r_ij = x[i] - x[j]
                distance = np.linalg.norm(r_ij)

                if distance != 0:
                    # Gravitational acceleration from body j on body i
                    a_i += -G * m[j] * r_ij / distance ** 3

        a[i] = a_i

    return a


def euler_method(m, W, h):
    """ Perform a single Euler step to update positions (x) and velocities (v). """

    # Split W into positions x and velocities v
    x = W[:, :3]    # First 3 columns are positions
    v = W[:, 3:]    # Next 3 columns are velocities

    # Compute accelerations based on current positions
    a = find_acceleration(m, x)

    # Update positions and velocities using the Euler method
    x_new = x + h * v   # Update positions
    v_new = v + h * a   # Update velocities

    # Combine new positions and velocities into W_new
    W_new = np.hstack((x_new, v_new))

    return W_new


# Simulate over time
num_steps = int(1 * 365 * 24 * 60 * 60 / h_val)     # Running for 1 Earth-year
W_history = [W_vals.copy()]  # To store the state at each time step

for step in range(num_steps):
    W_vals = euler_method(mass, W_vals, h_val)
    W_history.append(W_vals.copy())

# Extract positions history for plotting
positions_history = np.array([state[:, :3] for state in W_history])


def plot_trajectories(p_h):
    """ Plot the trajectories of all particles in 3D space. """

    labels = ["Sun", "Earth", "Mars"]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the trajectory of each particle
    for i in range(N):
        x_traj = p_h[:, i, 0]/1.496e13   # X trajectory of particle i
        y_traj = p_h[:, i, 1]/1.496e13   # Y trajectory of particle i
        z_traj = p_h[:, i, 2]/1.496e13   # Z trajectory of particle i
        ax.plot(x_traj, y_traj, z_traj, label=labels[i])

    ax.set_xlabel('X Position / AU')
    ax.set_ylabel('Y Position / AU')
    ax.set_zlabel('Z Position / AU')
    ax.set_title('Orbital Trajectories of Earth and Mars')
    ax.legend()

    plt.show()


# Plot the trajectories of the particles
plot_trajectories(positions_history)

# Example output of final positions after the simulation
x_final = W_vals[:, :3]
v_final = W_vals[:, 3:]

print("Final positions after simulation:")
print(x_final)
