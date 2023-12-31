
# import numpy as py

# import settings

# CO + 1/2 O2 --> CO2
# C3H6 + 9/2 O2 --> 3 CO2 + 3 H2O
# NO + CO --> CO2 + 1/2 N2
# H2 + 1/2 O2 --> H2O

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode


def twc(t, u, alpha, dx, dirichlet_bc, neumann_bc):

    du2dt2 = np.zeros_like(u)

    # Interior points
    du2dt2[1:-1] = alpha * (u[:-2] - 2 * u[1:-1] + u[2:]) / dx ** 2
    # Neumann boundary condition
    du2dt2[-1] = (neumann_bc - ((u[-1] - u[-2]) / dx)) / dx

    return du2dt2


# def jacobian(t, T, alpha, dx, dirichlet_bc, neumann_bc):
#
#     J = np.zeros((num_points, num_points))
#
#     # Interior points
#     J[1:-1, 1:-1] = -2 * alpha / dx ** 2
#
#     # Diagonal elements
#     J[1:-1, 1:-1] -= np.roll(J[1:-1, 1:-1], 1, axis=1)
#     J[1:-1, 1:-1] -= np.roll(J[1:-1, 1:-1], -1, axis=1)
#
#     # Neumann boundary condition
#     J[-1, -1] = -1 / dx
#
#     return J


def solve_twc(length, num_points, time_steps, alpha, dirichlet_bc, neumann_bc):
    dx = length / (num_points - 1)
    x = np.linspace(0, length, num_points)

    # Initial condition
    T0 = np.full(num_points, 0)
    T0[0] = dirichlet_bc

    # Set up ODE solver with Jacobian
    # r = ode(transient_heated_rod, jacobian).set_integrator('vode', method='bdf')
    r = ode(twc).set_integrator('vode', method='bdf')
    r.set_initial_value(T0)
    r.set_f_params(alpha, dx, dirichlet_bc, neumann_bc)
    # r.set_jac_params(alpha, dx, dirichlet_bc, neumann_bc)

    # Time-stepping loop
    results = [T0]
    for _ in range(time_steps):
        r.integrate(r.t + 1)  # Integrate to the next time step
        results.append(r.y)

    # Visualize results
    results = np.array(results).T
    for i in range(0, time_steps + 1, max(1, time_steps // 10)):
        plt.plot(x, results[:, i], label=f'Time Step {i}')

    plt.xlabel('Position (m)')
    plt.ylabel('Temperature')
    plt.legend()
    plt.show()


# Parameters
length = 1.0  # Length of the rod
num_points = 100  # Number of spatial points
time_steps = 1000  # Number of time steps
alpha = 0.01  # Thermal diffusivity
dirichlet_bc = 100.0  # Dirichlet boundary condition at x=0
neumann_bc = 100.0  # Neumann boundary condition at x=length

# Solve the transient heated rod problem using SciPy ode solver with Jacobian
solve_twc(length, num_points, time_steps, alpha, dirichlet_bc, neumann_bc)

