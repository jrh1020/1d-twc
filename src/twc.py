
import settings

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode


def twc(t, y, alpha, dx, T_inlet, wall_flux, num_z, num_y):

    print("Executing t = %.4f" % (t))

    #
    # VARIABLE RECALL AND PREP
    #

    # Temperature
    Ts     = y[:num_z]
    dTsdt = np.zeros_like(Ts)

    Tg     = y[num_z:2 * num_z]
    dTgdt = np.zeros_like(Tg)

    # Gas Species
    xCOg   = y[2 * num_z:3 * num_z]
    xO2g   = y[3 * num_z:4 * num_z]
    xCO2g  = y[4 * num_z:5 * num_z]
    xC3H6g = y[5 * num_z:6 * num_z]
    xH2Og  = y[6 * num_z:7 * num_z]
    xNOg   = y[7 * num_z:8 * num_z]
    xH2g   = y[8 * num_z:9 * num_z]
    xN2g   = y[9 * num_z:10 * num_z]

    dxCOgdt   = np.zeros_like(xCOg)
    dxO2gdt   = np.zeros_like(xO2g)
    dxCO2gdt  = np.zeros_like(xCO2g)
    dxC3H6gdt = np.zeros_like(xC3H6g)
    dxH2Ogdt  = np.zeros_like(xH2Og)
    dxNOgdt   = np.zeros_like(xNOg)
    dxH2gdt   = np.zeros_like(xH2g)
    dxN2gdt   = np.zeros_like(xN2g)

    # Solid-phase (washcoat) Species
    xCOw   = np.reshape(y[10 * num_z:(10 + num_y) * num_z], (num_z, num_y))
    xO2w   = np.reshape(y[(10+num_y) * num_z:(10 + 2*num_y) * num_z], (num_z, num_y))
    xCO2w  = np.reshape(y[(10+2*num_y) * num_z:(10 + 3*num_y) * num_z], (num_z, num_y))
    xC3H6w = np.reshape(y[(10+3*num_y) * num_z:(10 + 4*num_y) * num_z], (num_z, num_y))
    xH2Ow  = np.reshape(y[(10+4*num_y) * num_z:(10 + 5*num_y) * num_z], (num_z, num_y))
    xNOw   = np.reshape(y[(10+5*num_y) * num_z:(10 + 6*num_y) * num_z], (num_z, num_y))
    xH2w   = np.reshape(y[(10+6*num_y) * num_z:(10 + 7*num_y) * num_z], (num_z, num_y))
    xN2w   = np.reshape(y[(10+7*num_y) * num_z:(10 + 8*num_y) * num_z], (num_z, num_y))

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
    dTsdt[0]    = alpha * (wall_flux - ((Ts[0] - Ts[1]) / dx)) / dx
    dTsdt[1:-1] = alpha * (Ts[:-2] - 2 * Ts[1:-1] + Ts[2:]) / dx ** 2
    dTsdt[-1]   = alpha * (wall_flux - ((Ts[-1] - Ts[-2]) / dx)) / dx

    # Gas-phase temperature


    #
    # RETURN
    #
    return np.concatenate((dTsdt,
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
            dxN2wdt.flatten()))


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


def solve_twc(length, num_z, num_y, time_steps, alpha, inlet_temp, boundary_wall_flux, y0):
    dz = length / (num_z - 1)
    z = np.linspace(0, length, num_z)

    # Set up ODE solver with Jacobian
    r = (ode(lambda t, y: twc(t, y, alpha, dz, inlet_temp, boundary_wall_flux, num_z, num_y))
         .set_integrator('vode', method='bdf'))
    r.set_initial_value(y0)

    # Time-stepping loop
    y = [y0]
    t = [0]
    for _ in range(time_steps):
        r.integrate(r.t + 1)  # Integrate to the next time step
        y.append(r.y)
        t.append(r.t)

    # Visualize results
    # plt.plot(0.0, y0[:num_z])
    for i in range(len(t)):  # range(0, time_steps + 1, max(1, time_steps // 10)):
        plt.plot(z, y[i][:num_z], label=('t = %.2fs' % (t[i])))

    plt.xlabel('Position (m)')
    plt.ylabel('Temperature (K)')
    plt.legend()
    plt.show()


if __name__ == '__main__':

    num_z = settings.num_z
    num_y = settings.num_y
    y_initial = np.concatenate((
                                np.full((num_z, 1), settings.T0),  # Ts
                                np.full((num_z, 1), settings.T0),  # Tg

                                np.full((num_z, 1), settings.c0[0]),  # CO g
                                np.full((num_z, 1), settings.c0[1]),  # O2 g
                                np.full((num_z, 1), settings.c0[2]),  # CO2 g
                                np.full((num_z, 1), settings.c0[3]),  # C3H6 g
                                np.full((num_z, 1), settings.c0[4]),  # H2O g
                                np.full((num_z, 1), settings.c0[5]),  # NO g
                                np.full((num_z, 1), settings.c0[6]),  # H2 g
                                np.full((num_z, 1), settings.c0[7]),  # N2 g

                                np.full((num_z*num_y, 1), settings.c0[0]),  # CO wc
                                np.full((num_z*num_y, 1), settings.c0[1]),  # O2 wc
                                np.full((num_z*num_y, 1), settings.c0[2]),  # CO2 wc
                                np.full((num_z*num_y, 1), settings.c0[3]),  # C3H6 wc
                                np.full((num_z*num_y, 1), settings.c0[4]),  # H2O wc
                                np.full((num_z*num_y, 1), settings.c0[5]),  # NO wc
                                np.full((num_z*num_y, 1), settings.c0[6]),  # H2 wc
                                np.full((num_z*num_y, 1), settings.c0[7])  # N2 wc
                              ))

    # Solve the transient heated rod problem using SciPy ode solver
    solve_twc(settings.length,
              settings.num_z,
              settings.num_y,
              settings.time_steps,
              settings.alpha,
              settings.inlet_temp,
              settings.boundary_wall_flux,
              y_initial)
