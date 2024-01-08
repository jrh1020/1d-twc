
import settings

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode


def twc(t, y, alpha, dx, T_inlet, wall_flux, num_points, num_y):

    #
    # VARIABLE RECALL AND PREP
    #

    # Temperature
    Ts     = y[:num_points]
    dTsdt = np.zeros_like(Ts)

    Tg     = y[num_points + 1:2 * num_points]
    dTgdt = np.zeros_like(Tg)

    # Gas Species
    xCOg   = y[2*num_points+1:3*num_points]
    xO2g   = y[3*num_points+1:4*num_points]
    xCO2g  = y[4*num_points+1:5*num_points]
    xC3H6g = y[5*num_points+1:6*num_points]
    xH2Og  = y[6*num_points+1:7*num_points]
    xNOg   = y[7*num_points+1:8*num_points]
    xH2g   = y[8*num_points+1:9*num_points]
    xN2g   = y[9*num_points+1:10*num_points]

    dxCOgdt   = np.zeros_like(xCOg)
    dxO2gdt   = np.zeros_like(xO2g)
    dxCO2gdt  = np.zeros_like(xCO2g)
    dxC3H6gdt = np.zeros_like(xC3H6g)
    dxH2Ogdt  = np.zeros_like(xH2Og)
    dxNOgdt   = np.zeros_like(xNOg)
    dxH2gdt   = np.zeros_like(xH2g)
    dxN2gdt   = np.zeros_like(xN2g)

    # Solid-phase (washcoat) Species
    xCOw   = np.reshape(y[10*num_points+1:(10+num_y)*num_points], (num_points, num_y))
    xO2w   = np.reshape(y[(10+num_y)*num_points+1:(10+2*num_y)*num_points], (num_points, num_y))
    xCO2w  = np.reshape(y[(10+2*num_y)*num_points+1:(10+3*num_y)*num_points], (num_points, num_y))
    xC3H6w = np.reshape(y[(10+3*num_y)*num_points+1:(10+4*num_y)*num_points], (num_points, num_y))
    xH2Ow  = np.reshape(y[(10+4*num_y)*num_points+1:(10+5*num_y)*num_points], (num_points, num_y))
    xNOw   = np.reshape(y[(10+5*num_y)*num_points+1:(10+6*num_y)*num_points], (num_points, num_y))
    xH2w   = np.reshape(y[(10+6*num_y)*num_points+1:(10+7*num_y)*num_points], (num_points, num_y))
    xN2w   = np.reshape(y[(10+7*num_y)*num_points+1:(10+8*num_y)*num_points], (num_points, num_y))

    dxCOwdt   = np.zeros_like(xCOw)
    dxO2wdt   = np.zeros_like(xO2w)
    dxCO2wdt  = np.zeros_like(xCO2w)
    dxC3H6wdt = np.zeros_like(xC3H6w)
    dxH2Owdt  = np.zeros_like(xH2Ow)
    dxNOwdt   = np.zeros_like(xNOw)
    dxH2wdt   = np.zeros_like(xH2w)
    dxN2wdt   = np.zeros_like(xN2w)

    #
    # CALCULATIONS
    #

    # Solid-phase temperature

    # Neumann boundary condition (ghost point method)
    dTsdt[0] = (wall_flux - ((Ts[0] - Ts[1]) / dx)) / dx
    dTsdt[1:-1] = alpha * (Ts[:-2] - 2 * Ts[1:-1] + Ts[2:]) / dx ** 2
    dTsdt[-1] = (wall_flux - ((Ts[-1] - Ts[-2]) / dx)) / dx

    #
    # RETURN
    #
    return [dTsdt,
            dTgdt,
            dxCOgdt,
            dxO2gdt,
            dxCO2gdt,
            dxC3H6gdt,
            dxH2Ogdt,
            dxNOgdt,
            dxH2gdt,
            dxN2gdt,
            dxCOwdt.flatten(),
            dxO2wdt.flatten(),
            dxCO2wdt.flatten(),
            dxC3H6wdt.flatten(),
            dxH2Owdt.flatten(),
            dxNOwdt.flatten(),
            dxH2wdt.flatten(),
            dxN2wdt.flatten()]


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


def solve_twc(length, num_points, time_steps, alpha, dirichlet_bc, neumann_bc, y_initial):
    dx = length / (num_points - 1)
    x = np.linspace(0, length, num_points)

    # Initial condition
    y0 = np.full(num_points, y_initial)

    # Set up ODE solver with Jacobian
    # r = ode(transient_heated_rod, jacobian).set_integrator('vode', method='bdf')
    r = ode(twc).set_integrator('vode', method='bdf')
    r.set_initial_value(y0)
    r.set_f_params(alpha, dx, dirichlet_bc, neumann_bc, num_points)
    # r.set_jac_params(alpha, dx, dirichlet_bc, neumann_bc)

    # Time-stepping loop
    results = [y0]
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


# Solve the transient heated rod problem using SciPy ode solver
solve_twc(settings.length,
          settings.num_points,
          settings.time_steps,
          settings.alpha,
          settings.dirichlet_bc,
          settings.neumann_bc,
          [settings.T0,
                  settings.T0])
