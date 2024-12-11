import numpy as np
import matplotlib.pyplot as plt


# Define parameters
nx = 2000
nsteps = 0
cartesian = False
pi = np.pi
gamma, t, dt, dt12, tmax, q, fmratio, cflfactor, dfactor, efac, il, istep = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

# Define arrays
fm = np.zeros(nx)
dm = np.zeros(nx)
dm12 = np.zeros(nx)
r = np.zeros(nx)
dr = np.zeros(nx)
dr12 = np.zeros(nx)
rold = np.zeros(nx)
v = np.zeros(nx)
a = np.zeros(nx)
at12 = np.zeros(nx)
ak12 = np.zeros(nx)
u = np.zeros(nx)
p = np.zeros(nx)
rho = np.zeros(nx)
eps = np.zeros(nx)
w = np.zeros(nx)
wt12 = np.zeros(nx)
aux = np.zeros(nx)


def read_input(filename, case_number=1):    # Helper function for reading input parameters
    # Reads input parameters from file and assigns them to common_vars based on the chosen case
    global nsteps, tmax, efac, il, q, fmratio

    # Each case in the file has 8 lines (1 description + 6 values + 1 blank line)
    start_line = (case_number - 1) * 8 + 1
    with open(filename, 'r') as file:
        lines = file.readlines()
        # Extract the specific case's parameters
        params = [float(lines[start_line + i].split(',')[0].strip()) for i in range(6)]

    # Assign values to the common_vars dictionary
    nsteps = int(params[0])
    tmax = int(params[1])
    efac = params[2]
    il = int(params[3])
    q = params[4]
    fmratio = params[5]


def inicond():
    global cflfactor, dfactor, gamma, il, fmratio, pi, dm12, rho, eps, efac, p, u, fm, dm, r, v, cartesian, a, at12, \
        dr12, t
    # Initialize the physical conditions based on initial parameters
    cflfactor = 0.1
    dfactor = 1.5 ** 2
    gamma = 5.0 / 3.0

    # High pressure ejecta
    rhoej = ((float(nx) / float(il)) ** 3 - 1.0) * fmratio
    dmej = 4.0 / 3.0 * pi * (float(il) / float(nx)) ** 3 * rhoej / float(il)
    for ix in range(il):
        dm12[ix] = dmej
        rho[ix] = rhoej
        eps[ix] = efac
        p[ix] = (gamma - 1.0) * rho[ix] * eps[ix]
        u[ix] = 0.0

    # Low pressure ambient medium
    rhoamb = 1.0
    dmamb = 4.0 / 3.0 * pi * (1.0 - (float(il) / float(nx)) ** 3) / float(nx - il) * rhoamb
    for ix in range(il, nx):
        dm12[ix] = dmamb
        rho[ix] = rhoamb
        eps[ix] = 1.0
        p[ix] = (gamma - 1.0) * rho[ix] * eps[ix]
        u[ix] = 0.0

    # Mass and volume initialisation
    fm[0] = dm12[0]
    for ix in range(1, nx):
        fm[ix] = fm[ix - 1] + dm12[ix]
        dm[ix] = 0.5 * (dm12[ix] + dm12[ix - 1])

    # Radius and area setup
    r[0] = 0.0
    v[0] = 0.0
    for ix in range(1, nx):
        v[ix] = v[ix - 1] + dm12[ix - 1] / rho[ix - 1]
        if cartesian:
            r[ix] = v[ix]
            a[ix] = 1.0
            at12[ix] = 1.0
        else:
            r[ix] = (v[ix] / (4.0 / 3.0 * pi)) ** (1.0 / 3.0)
            a[ix] = 4.0 * pi * r[ix] ** 2
    for ix in range(nx - 1):
        dr12[ix] = r[ix + 1] - r[ix]

    t = 0.0


def tmsc():
    global dr12, u, gamma, eps, cflfactor, t, tmax, at12, v, dfactor, dt12, dt

    dtc = 1.0e30
    for ix in range(nx - 1):
        dtc = min(dtc, dr12[ix] / (abs(u[ix]) + np.sqrt(gamma * eps[ix])))
    dtc *= cflfactor
    if t + dtc > tmax:
        dtc = tmax - t

    # Diffusion limit
    dtd = 1.0e-30
    for ix in range(nx - 1):
        dtd = max(dtd, float(abs(at12[ix + 1] * u[ix + 1] - at12[ix] * u[ix]) / (v[ix + 1] - v[ix])))
    dtd = 0.5 / (dtd * dfactor)

    dt12 = min(dtc, dtd)
    dt = 0.5 * (dt12 + dtc)


def run_simulation(case_number):
    # Read input parameters and initialize conditions
    read_input("lh1.txt", case_number)
    inicond()
    global istep, nsteps, u, a, p, dt, dm, w, ak12, rold, r, dt12, dr12, v, at12, pi, rho, dm12, q, aux, eps, gamma, \
        t, tmax

    # Open output files
    with open("lh1.out", "w") as out1, open("lh2.out", "w") as out2:
        etot0 = 1.0  # Total initial energy to normalize future energies

        # Loop over time steps
        for istep in range(1, nsteps + 1):
            tmsc()  # Call function to determine time step size `dt`

            # Update velocities
            for ix in range(1, nx):
                u[ix] = (u[ix] - a[ix] * (p[ix] - p[ix - 1]) * dt / dm[ix]
                         - 0.5 * (w[ix] * (3.0 * ak12[ix] - a[ix]) - w[ix - 1] *
                         (3.0 * ak12[ix - 1] - a[ix])) * dt / dm[ix])

            # Apply boundary conditions
            if cartesian:
                u[0] = u[1]
            else:
                u[0] = 0.0

            # Update radii, surfaces, and volumes
            rold[:] = r
            r[:] = rold + u * dt12
            dr12[: nx - 1] = r[1:] - r[:-1]

            if cartesian:
                v[:] = r
            else:
                at12[:] = 4.0 * pi * (0.5 * (r + rold)) ** 2
                a[:] = 4.0 * pi * r ** 2
                v[:] = 4.0 / 3.0 * pi * r ** 3
                ak12[: nx - 1] = 0.5 * (at12[1:] + at12[:-1])

            # Update densities
            rho[: nx - 1] = dm12[:-1] / (v[1:] - v[:-1])
            rho[-1] = rho[-2]

            # Apply artificial viscosity
            for ix in range(nx - 1):
                if u[ix + 1] > u[ix]:
                    w[ix] = 0.0
                else:
                    w[ix] = - q ** 2 * rho[ix] * abs(u[ix + 1] - u[ix]) * (
                            u[ix + 1] * (1.0 - at12[ix + 1] / (3.0 * ak12[ix])) - u[ix] *
                            (1.0 - at12[ix] / (3.0 * ak12[ix])))

            # Update internal energies and pressures
            for ix in range(nx - 1):
                aux[ix] = eps[ix] - p[ix] * (at12[ix + 1] * u[ix + 1] - at12[ix] * u[ix]) * dt12 / dm12[ix]
                p[ix] = 0.5 * (p[ix] + (gamma - 1.0) * rho[ix] * aux[ix])
                eps[ix] = eps[ix] - p[ix] * (at12[ix + 1] * u[ix + 1] - at12[ix] * u[ix]) * dt12 / dm12[ix]

            # Contribution from artificial viscosity
            for ix in range(nx - 1):
                eps[ix] -= 0.5 * w[ix] * dt12 / dm12[ix] * (
                    u[ix + 1] * (3.0 * ak12[ix] - at12[ix + 1]) - u[ix] * (3.0 * ak12[ix] - at12[ix]))

            # Update pressures for remaining shells
            p[: nx - 1] = (gamma - 1.0) * rho[: nx - 1] * eps[: nx - 1]
            p[-1] = p[-2]
            eps[-1] = eps[-2]

            # Update time
            t += dt12

            # Energy conservation check
            ethe = np.sum(eps[1:nx - 1] * dm[1:nx - 1])
            ekin = 0.5 * np.sum((0.5 * (u[2:] + u[1:nx - 1])) ** 2 * dm[1:nx - 1])
            etot = ethe + ekin
            if istep == 1:
                etot0 = etot
            etot /= etot0
            ethe /= etot0
            ekin /= etot0

            # Write energies every 10 steps to out2 and every 100 steps to console
            if istep % 10 == 0:
                out2.write(f"{istep:5d} {t:10.3e} {etot:10.3e} {ethe:10.3e} {ekin:10.3e}\n")
            if istep % 100 == 0:
                print(f"{istep:5d}; t: {t:10.3e}; etot: {etot:10.3e}; eth: {ethe:10.3e}; ekin: {ekin:10.3e}")

            # End time step check
            if t >= tmax:
                break

        # Final printout
        umax = max(u)
        rhomax = max(rho)
        pmax = max(p)
        emax = max(eps)
        wmax = max(w)

        # Normalize final values and write to file
        epsi = 1.0e-20
        u /= (epsi + umax)
        rho /= (epsi + rhomax)
        p /= (epsi + pmax)
        eps /= (epsi + emax)
        w /= (epsi + wmax)

        for ix in range(nx):
            out1.write(
                f"{ix + 1:5d} {fm[ix]:12.4e} {r[ix] / r[nx - 1]:12.4e} {u[ix]:12.4e} "
                f"{rho[ix]:12.4e} {p[ix]:12.4e} {eps[ix]:12.4e} {w[ix]:12.4e}\n")


def plot_results():
    """Plots the simulation results from output files."""

    # Load final model data from 'lh1.out'
    data = np.loadtxt("lh1.out")
    indices_plot = data[:, 0]           # Radial index
    fm_plot = data[:, 1]                # Mass fraction
    r_norm_plot = data[:, 2]            # Normalized radius
    u_plot = data[:, 3]                 # Velocity
    rho_plot = data[:, 4]               # Density
    p_plot = data[:, 5]                 # Pressure
    eps_plot = data[:, 6]               # Internal energy
    w_plot = data[:, 7]                 # Artificial viscosity

    # Load energies over time from 'lh2.out'
    time_data = np.loadtxt("lh2.out")
    time_steps_plot = time_data[:, 0]   # Time step index
    time_plot = time_data[:, 1]         # Time
    etot_plot = time_data[:, 2]         # Total energy
    ethe_plot = time_data[:, 3]         # Thermal energy
    ekin_plot = time_data[:, 4]         # Kinetic energy

    # Create subplots for each quantity
    fig, axs = plt.subplots(4, 2, figsize=(14, 20))
    fig.suptitle("Simulation Results")

    # Plot velocity
    axs[0, 0].plot(r_norm_plot, u_plot, label="Velocity (u)", color="b")
    axs[0, 0].set_xlabel("Normalized Radius")
    axs[0, 0].set_ylabel("Velocity (u)")
    axs[0, 0].set_title("Velocity Profile")
    axs[0, 0].legend()

    # Plot density
    axs[0, 1].plot(r_norm_plot, rho_plot, label="Density (rho)", color="g")
    axs[0, 1].set_xlabel("Normalized Radius")
    axs[0, 1].set_ylabel("Density (rho)")
    axs[0, 1].set_title("Density Profile")
    axs[0, 1].legend()

    # Plot pressure
    axs[1, 0].plot(r_norm_plot, p_plot, label="Pressure (p)", color="r")
    axs[1, 0].set_xlabel("Normalized Radius")
    axs[1, 0].set_ylabel("Pressure (p)")
    axs[1, 0].set_title("Pressure Profile")
    axs[1, 0].legend()

    # Plot internal energy
    axs[1, 1].plot(r_norm_plot, eps_plot, label="Internal Energy (eps)", color="purple")
    axs[1, 1].set_xlabel("Normalized Radius")
    axs[1, 1].set_ylabel("Internal Energy (eps)")
    axs[1, 1].set_title("Internal Energy Profile")
    axs[1, 1].legend()

    # Plot artificial viscosity
    axs[2, 0].plot(r_norm_plot, w_plot, label="Artificial Viscosity (w)", color="orange")
    axs[2, 0].set_xlabel("Normalized Radius")
    axs[2, 0].set_ylabel("Artificial Viscosity (w)")
    axs[2, 0].set_title("Artificial Viscosity Profile")
    axs[2, 0].legend()

    # Plot total energy over time
    axs[2, 1].plot(time_plot, etot_plot, label="Total Energy (etot)", color="black")
    axs[2, 1].set_xlabel("Time")
    axs[2, 1].set_ylabel("Total Energy")
    axs[2, 1].set_title("Total Energy Over Time")
    axs[2, 1].legend()

    # Plot thermal energy over time
    axs[3, 0].plot(time_plot, ethe_plot, label="Thermal Energy (ethe)", color="brown")
    axs[3, 0].set_xlabel("Time")
    axs[3, 0].set_ylabel("Thermal Energy")
    axs[3, 0].set_title("Thermal Energy Over Time")
    axs[3, 0].legend()

    # Plot kinetic energy over time
    axs[3, 1].plot(time_plot, ekin_plot, label="Kinetic Energy (ekin)", color="cyan")
    axs[3, 1].set_xlabel("Time")
    axs[3, 1].set_ylabel("Kinetic Energy")
    axs[3, 1].set_title("Kinetic Energy Over Time")
    axs[3, 1].legend()

    # Adjust layout and spacing for the main figure
    plt.tight_layout(rect=(0, 0.03, 1, 0.95))

    # Create a separate figure for the combined plot of u, rho, p, and eps
    fig_combined, ax_combined = plt.subplots(figsize=(10, 6))
    ax_combined.plot(r_norm_plot, u_plot, label="Velocity (u)", color="b", linestyle='-')
    ax_combined.plot(r_norm_plot, rho_plot, label="Density (rho)", color="g", linestyle='--')
    ax_combined.plot(r_norm_plot, p_plot, label="Pressure (p)", color="r", linestyle='-.')
    ax_combined.plot(r_norm_plot, w_plot, label="Artificial Viscosity (w)", color="purple", linestyle=':')
    ax_combined.set_xlabel("Normalized Radius")
    ax_combined.set_ylabel("Magnitude")
    ax_combined.set_title("Combined Profile (Velocity, Density, Pressure, Artificial Viscosity)")
    ax_combined.legend()

    # Show all figures
    plt.show()


run_simulation(1)
plot_results()
