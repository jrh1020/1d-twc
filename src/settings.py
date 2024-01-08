
# Species:
# 1: CO
# 2: O2
# 3: CO2
# 4: C3H6
# 5: H2O
# 6: NO
# 7: H2
# 8: N2

# CO + 1/2 O2 --> CO2
# C3H6 + 9/2 O2 --> 3 CO2 + 3 H2O
# NO + CO --> CO2 + 1/2 N2
# H2 + 1/2 O2 --> H2O

# Parameters
length = 0.1  # Length of the rod (m)
num_points = 100  # Number of spatial points
num_y = 5  # Number of washcoat discretizations
time_steps = 1000  # Number of time steps
rho = 8850  # Density (kg/m3)
lam = 398  # Thermal conductivity (W/m/K)
cp = 390  # Heat capacity (J/kg/K)
alpha = lam / rho / cp
# alpha = 1.11e-4  # Thermal diffusivity
dirichlet_bc = 373.15  # Dirichlet boundary condition at x=0
neumann_bc = 100.0  # Neumann boundary condition at x=length

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