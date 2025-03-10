# Import Modules
import modules

# Code for Turbojet Cycle Analysis
# Functions

def compressor_inlet(gamma, mach, altitude):
    # Get atmospheric conditions at given altitude in m
    atmos_con = modules.atmos(altitude)
    
    # calculate compressor inlet temperature
    T_02 = atmos_con.T[0] * (1 + ((gamma - 1)/2 * mach**2))

    # assumption, want this to be calculated eventually, maybe can do iteration similar to comp_outlet
    diffuser_eff = 0.97

    # calculate compressor inlet pressure
    P_02 = atmos_con.P[0] * (1 + (diffuser_eff * (T_02/atmos_con.T[0] - 1))) ** (gamma/(gamma - 1))

    # Compressor Inlet Temperature, Pressure and Atmospheric conditions
    return atmos_con.T[0], atmos_con.P[0], T_02, P_02, atmos_con

def compressor_outlet(pressure_ratio, P_02, T_02, atmos_con, gamma):
    # obtain compressor outlet pressure
    P_03 = pressure_ratio * P_02

    # adiabatic compressor efficiency
    adiabatic_comp_eff = ((P_02/atmos_con.P[0])**((gamma-1)/gamma)-1)/((T_02/atmos_con.T[0])-1)

    # specific heat ratio for the compression process
    max_iter = 50
    tol = 1e-5

    gamma_c = gamma # initial guess based on air = 1.4

    # obtain compressor outlet temperature through iteration

    for i in range(max_iter):
        # find compressor outlet temperature
        T_03 = T_02 * (1 + 1/(adiabatic_comp_eff) * (pressure_ratio**((gamma_c - 1)/gamma_c) - 1))

        T_avg = (T_02 + T_03)/2
        P_avg = (P_02 + P_03)/2

        # use CoolProps to find the specific heat and volume capacity of air given conditions
        Cp = modules.cp.PropsSI('Cpmass', 'T', T_avg, 'P', P_avg, 'Air')
        Cv = modules.cp.PropsSI('Cvmass', 'T', T_avg, 'P', P_avg, 'Air')
        gamma_c_new = Cp/Cv

        # check tolerance of 
        if abs(gamma_c_new - gamma_c) < tol:
            gamma_c = gamma_c_new
            break
        gamma_c = gamma_c_new

    # Compressor outlet Temperature, Pressure, Specific Heat Constant Pressure
    return T_03, P_03, Cp

def fuel_to_air_ratio(T_04, T_03, Q_R, P_03, Cp):
    # find fuel to air ratio using given conditions
    fuel_air_ratio = ((T_04/T_03) - 1) / ((Q_R/(Cp * T_03)) - (T_04/T_03))

    # fuel to air combustion ratio
    return fuel_air_ratio

def turbine_outlet(T_04, T_03, T_02, P_04): 
    # assume turbine adaibatic efficiency
    turbine_adiabatic_eff = 0.9

    # Simplified method to find Turbine Outlet Temperature
    T_05 = T_04 - (T_03 - T_02)

    # specific heat ratio for the turbine process
    max_iter = 50
    tol = 1e-5

    gamma_t = gamma # initial guess based on air = 1.4

    # obtain turbine outlet pressure through iteration

    for i in range(max_iter):
        # find turbine outlet pressure
        P_05 = P_04 * (1 - 1/(turbine_adiabatic_eff) * (1 - (T_05/T_04))) ** (gamma_t/(gamma_t - 1))

        T_avg = (T_05 + T_04)/2
        P_avg = (P_05 + P_04)/2

        # use CoolProps to find the specific heat and volume capacity of air given conditions
        Cp = modules.cp.PropsSI('Cpmass', 'T', T_avg, 'P', P_avg, 'Air')
        Cv = modules.cp.PropsSI('Cvmass', 'T', T_avg, 'P', P_avg, 'Air')
        gamma_t_new = Cp/Cv

        # check tolerance of 
        if abs(gamma_t_new - gamma_t) < tol:
            gamma_t = gamma_t_new
            break
        gamma_t = gamma_t_new
    
    # Nozzle inlet Temperature and Pressure
    return T_05, P_05

def nozzle_velocity(T_06, P_06, R, atmos_con):
    # nozzle gamma
    gamma_n = 1.3
    
    nozzle_adiabatic_eff = 0.98
    # exhaust gas exit velocity
    exit_velocity = (2 * nozzle_adiabatic_eff * ((gamma_n/(gamma_n - 1))) * R * T_06 * (1 - (atmos_con.P[0] / P_06) ** ((gamma_n - 1)/gamma_n))) ** 0.5

    # gas exhaust exit velocity
    return exit_velocity

def fuel_flow(T_05, P_05, mach, T_01, P_01, R, gamma, exit_velocity, area_inlet, area_exit, fuel_air_ratio):
    # calculate densities from atmospheric conditions
    rho = modules.cp.PropsSI('D', 'T', T_01, 'P', P_01, 'Air')
    rho_e = modules.cp.PropsSI('D', 'T', T_05, 'P', P_05, 'Air')

    # calculate vehicle speed in terms of m/s
    a = (gamma * R * T_01) ** 0.5
    inlet_velocity = a * mach

    # calculate fuel flow rate in kg/s
    mass_fuel_flow_rate = rho_e * exit_velocity * area_exit - rho * inlet_velocity * area_inlet
    mass_air_flow_rate = rho * inlet_velocity * area_inlet
    mass_exit_flow_rate = mass_fuel_flow_rate + mass_air_flow_rate

    return mass_fuel_flow_rate, mass_air_flow_rate, mass_exit_flow_rate, inlet_velocity

def thrust(mass_air_flow_rate, fuel_air_ratio, exit_velocity, inlet_velocity, P_05, P_01, area_exit):
    # calculate thrust using conservation of momentum
    thrust = mass_air_flow_rate*((1 + fuel_air_ratio) * exit_velocity - inlet_velocity) + (P_05 - P_01) * area_exit

    return thrust


# Modelled after the -, for now, want to find some optimisation formula
engine_TJ_mass = 1  # kg, dry mass, could become user defined input
max_TJ_thrust = 1 # Newtons, could become user defined input
max_TJ_mach = 2 # to be revised
gamma = 1.4 # for air, pre combustion
pressure_ratio = 10 # for now, want this to become an user defined input eventually
max_TJ_temp = 1700 # Kelvin
Q_R = 48E6 # J/kg
R = 287 # J / (kg * K)
area_inlet = 0.7    # m², typical turbojet inlet area
area_exit = 0.2     # m², typical exhaust/nozzle throat area
area_ratio = area_exit / area_inlet  # approximately 0.2857


mach = 1
altitude = 10

# inputs: gamma, mach number, altitude in km
T_01, P_01, T_02, P_02, atmos_con = compressor_inlet(gamma, mach, altitude)

# inputs: pressure ratio (variable), compressor inlet pressure, compressor inlet temperature, atmospheric conditions, gamma
T_03, P_03, Cp = compressor_outlet(pressure_ratio, P_02, T_02, atmos_con, gamma)

# inputs: maximum temperature, compressor outlet temperature, specific energy release of the fuel, compressor outlet pressure, specific heat
fuel_air_ratio = fuel_to_air_ratio(max_TJ_temp, T_03, Q_R, P_03, Cp)

# inputs: maximum temperature, compressor outlet temperature, compressor inlet temperature, compressor outlet pressure
T_05, P_05 = turbine_outlet(max_TJ_temp, T_03, T_02, P_03)

# inputs: nozzle exit temperature, nozzle exit pressure, specific gas constant for air
exit_velocity = nozzle_velocity(T_05, P_05, R, atmos_con)

# inputs: nozzle exit temperature, nozzle exit pressure, mach number, ambient temperature, ambient pressure, specific gas constant, gamma
#         nozzle exit velocity, inlet area, exhaust area, fuel to air ratio
mass_fuel_flow_rate, mass_air_flow_rate, mass_exit_flow_rate, inlet_velocity = fuel_flow(T_05, P_05, mach, T_01, P_01, R, gamma,
exit_velocity, area_inlet, area_exit, fuel_air_ratio)

# inputs: air mass flow rate, fuel to air ratio, nozzle exit velocity, inlet velocity, nozzle exit pressure, ambient pressure, nozzle area
thrust = thrust(mass_air_flow_rate, fuel_air_ratio, exit_velocity, inlet_velocity, P_05, P_01, area_exit)
print(f"The TJ produces {thrust:.2f}N of Thrust")

temperatures = [T_01, T_02, T_03, max_TJ_temp, T_05]
pressure = [P_01, P_02, P_03, P_03, P_05]

# Calculate specific volume and entropy
specific_volumes = [R * T / P for T, P in zip(temperatures, pressure)]
entropies = [modules.cp.PropsSI('S', 'T', T, 'P', P, 'Air') for T, P in zip(temperatures, pressure)]

# Plot results
modules.plt.figure()
modules.plt.plot(entropies, temperatures)
modules.plt.plot([entropies[0], entropies[-1]], [temperatures[0], temperatures[-1]],
                 ':', color='gray', linewidth=1, label='Start/End Connection')
modules.plt.xlabel('Entropy (J/kg·K)')
modules.plt.ylabel('Temperature (K)')
modules.plt.title('T-S Diagram')

modules.plt.figure()
modules.plt.plot(specific_volumes, pressure)
modules.plt.plot([specific_volumes[0], specific_volumes[-1]], [pressure[0], pressure[-1]],  # Connect start and end
                 ':', color='gray', linewidth=1, label='Start/End Connection')
modules.plt.xlabel('Specific Volume (m³/kg)')
modules.plt.ylabel('Pressure (Pa)')
modules.plt.title('P-v Diagram')

modules.plt.show()









