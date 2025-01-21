# Species:
# 1: CO
# 2: O2
# 3: CO2
# 4: C3H6
# 5: H2O
# 6: NO
# 7: H2
# 8: N2

c0 = [0.0,
      0.21,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.79]

# Diffusion coefficients (m2/s/K^1.75) from Fuller, Schettler, Giddings
D0 = [1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0]
Dexp = [1.75 for _ in range(len(D0))]  #,  # Fuller, Schettler, Giddings Exponents
        # 1.75,
        # 1.75,
        # 1.75,
        # 1.75,
        # 1.75,
        # 1.75,
        # 1.75]

# Reactions:
# 1: CO + 1/2 O2 --> CO2
# 2: C3H6 + 9/2 O2 --> 3 CO2 + 3 H2O
# 3: NO + CO --> CO2 + 1/2 N2
# 4: H2 + 1/2 O2 --> H2O

# Reaction coefficients
#       CO    O2    CO2  C3H6   H2O   NO    H2    N2
nu = [[-1.0, -0.5, +1.0, +0.0, +0.0, +0.0, +0.0, +0.0],  # 1
      [+0.0, -4.5, +3.0, -1.0, +3.0, +0.0, +0.0, +0.0],  # 2
      [-1.0, +0.0, +1.0, +0.0, +0.0, -1.0, +0.0, +0.5],  # 3
      [+0.0, -0.5, +0.0, +0.0, +1.0, +0.0, +0.0, +0.0]   # 4
      ]

# Arrhenius equation ccoefficients
Arr = [1.0,
       1.0,
       1.0,
       1.0]
Ea = [10.0,
      10.0,
      10.0,
      10.0]

# Heats of reactions (kJ/mol)
delH = [-283.0,
        -1926.3,
        -373.25,
        -241.8]

# Parameters
time_steps = 10  # Number of time steps

length = 0.1  # Length of the rod (m)
num_z = 100  # Number of spatial points
num_y = 5  # Number of washcoat discretizations

rho = 8850  # Density (kg/m3)
lam = 398  # Thermal conductivity (W/m/K)
cp = 390  # Heat capacity (J/kg/K)
alpha = lam / rho / cp

# alpha = 1.11e-4  # Thermal diffusivity
inlet_temp = 373.15  # Dirichlet boundary condition at x=0
boundary_wall_flux = 0.0  # Neumann boundary condition at x=length

# == Constants ==
Sh = 2.98
Nu = 2.98

# == Settings ==
# Geometry
cpsi = 900  # 1/in^2
t_wc = 25e-6  # m

# Boundary Conditions
v = 3  # m/s
Tgin = 600  # K

# Initial Conditions
T0 = 298.15  # K
