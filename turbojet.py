# ------------------ Import ------------------#
import modules

# ------------------ Functions ------------------#
def compressor_inlet(gamma, mach, altitude):
    """
    Calculate the conditions at the inlet of the compressor

    Parameters:
        gamma (float): Specific heat ratio of the air
        mach (float): Mach number
        altitude (float): Altitude in kilometers

    Returns:
        tuple: Contains:
            T_01 (float): Ambient temperature at the inlet
            P_01 (float): Ambient pressure at the inlet
            T_02 (float): Compressor inlet temperature
            P_02 (float): Compressor inlet pressure
    """
    atmos_con = modules.atmos(altitude)
    if mach == 0:
        T_01 = 288.15
        P_01 = atmos_con.P[0]
    else:
        # Get atmospheric conditions at given altitude in km
        T_01 = atmos_con.T[0]
        P_01 = atmos_con.P[0]
        
    # calculate compressor inlet temperature
    T_02 = T_01 * (1 + ((gamma - 1)/2 * mach ** 2))

    # assumption, want this to be calculated eventually, maybe can do iteration similar to comp_outlet
    diffuser_eff = 0.97

    # calculate compressor inlet pressure
    P_02 = P_01 * (1 + (diffuser_eff * (T_02 / T_01 - 1))) ** (gamma/(gamma - 1))

    # Compressor Inlet Temperature, Pressure and Atmospheric conditions
    return T_01, P_01, T_02, P_02

def compressor_outlet(pressure_ratio, P_02, T_02, gamma, T_01, P_01):
    """
    Calculate the conditions at the outlet of the compressor

    Parameters:
        pressure_ratio (float): Compressor pressure ratio
        P_02 (float): Compressor inlet pressure
        T_02 (float): Compressor inlet temperature
        gamma (float): Specific heat ratio of the air
        T_01 (float): Ambient temperature at the inlet
        P_01 (float): Ambient pressure at the inlet

    Returns:
        tuple: Contains:
            T_03 (float): Compressor outlet temperature
            P_03 (float): Compressor outlet pressure
            Cp (float): Specific heat at constant pressure
    """
    # obtain compressor outlet pressure
    P_03 = pressure_ratio * P_02

    if T_02 == T_01:
        adiabatic_comp_eff = 1
    else:
        # adiabatic compressor efficiency
        adiabatic_comp_eff = ((P_02 / P_01) ** ((gamma - 1) / gamma) - 1)/((T_02 / T_01 ) - 1)

    # specific heat ratio for the compression process
    max_iter = 50
    tol = 1e-5

    gamma_c = gamma  # initial guess based on air = 1.4

    # obtain compressor outlet temperature through iteration

    for i in range(max_iter):
        # find compressor outlet temperature
        T_03 = T_02 * (1 + 1 / (adiabatic_comp_eff) * (pressure_ratio **((gamma_c - 1)/ gamma_c) - 1))

        T_avg = (T_02 + T_03) / 2
        P_avg = (P_02 + P_03) / 2

        # use CoolProps to find the specific heat and volume capacity of air given conditions
        Cp = modules.cp.PropsSI('Cpmass', 'T', T_avg, 'P', P_avg, 'Air')
        Cv = modules.cp.PropsSI('Cvmass', 'T', T_avg, 'P', P_avg, 'Air')
        gamma_c_new = Cp / Cv

        # check tolerance of gamma_c to check for convergence
        if abs(gamma_c_new - gamma_c) < tol:
            gamma_c = gamma_c_new
            break
        gamma_c = gamma_c_new

    # Compressor outlet Temperature, Pressure, Specific Heat Constant Pressure
    return T_03, P_03, Cp

def fuel_to_air_ratio(T_04, T_03, Q_R, P_03, Cp):
    """
    Calculate the fuel-to-air ratio for the engine

    Parameters:
        T_04 (float): Maximum temperature (combustor exit temperature)
        T_03 (float): Compressor outlet temperature
        Q_R (float): Specific energy release of the fuel (J/kg)
        P_03 (float): Compressor outlet pressure
        Cp (float): Specific heat at constant pressure

    Returns:
        float: Fuel-to-air ratio
    """
    # find fuel to air ratio using given conditions, in
    fuel_air_ratio = ((T_04 / T_03) - 1) / ((Q_R / (Cp * T_03)) - (T_04 / T_03))

    # fuel to air combustion ratio
    return fuel_air_ratio

def turbine_outlet(T_04, T_03, T_02, P_04, gamma):
    """
    Calculate the turbine outlet conditions

    Parameters:
        T_04 (float): Maximum temperature (combustor exit temperature)
        T_03 (float): Compressor outlet temperature
        T_02 (float): Compressor inlet temperature
        P_04 (float): Compressor outlet pressure
        gamma (float): Specific heat ratio of the air

    Returns:
        tuple: Contains:
            T_05 (float): Turbine outlet temperature
            P_05 (float): Turbine outlet pressure
    """
    # assume turbine adaibatic efficiency
    turbine_adiabatic_eff = 0.9

    # Simplified method to find Turbine Outlet Temperature
    T_05 = T_04 - (T_03 - T_02)

    # specific heat ratio for the turbine process
    max_iter = 50
    tol = 1e-5

    gamma_t = gamma  # initial guess based on air = 1.4

    # obtain turbine outlet pressure through iteration
    for i in range(max_iter):
        # find turbine outlet pressure
        P_05 = P_04 * (1 - 1 / (turbine_adiabatic_eff) * (1 - (T_05 / T_04))) ** (gamma_t / (gamma_t - 1))

        T_avg = (T_05 + T_04)/2
        P_avg = (P_05 + P_04)/2

        # use CoolProps to find the specific heat and volume capacity of air given conditions
        Cp = modules.cp.PropsSI('Cpmass', 'T', T_avg, 'P', P_avg, 'Air')
        Cv = modules.cp.PropsSI('Cvmass', 'T', T_avg, 'P', P_avg, 'Air')
        gamma_t_new = Cp/Cv

        # check tolerance of gamma_t to check for convergence
        if abs(gamma_t_new - gamma_t) < tol:
            gamma_t = gamma_t_new
            break
        gamma_t = gamma_t_new
        
    # Nozzle inlet Temperature and Pressure
    return T_05, P_05

def nozzle_velocity(T_06, P_06, R, P_01):
    """
    Calculate the nozzle exit velocity

    Parameters:
        T_06 (float): Nozzle exit temperature
        P_06 (float): Nozzle exit pressure
        R (float): Specific gas constant for air
        P_01 (float): Ambient pressure at the inlet

    Returns:
        float: Nozzle exit velocity in m/s
    """
    # nozzle gamma
    gamma_n = 1.3

    nozzle_adiabatic_eff = 0.98
    # exhaust gas exit velocity
    exit_velocity = (2 * nozzle_adiabatic_eff * ((gamma_n / (gamma_n - 1))) * R * T_06 * (1 - (P_01 / P_06) ** ((gamma_n - 1)/ gamma_n))) ** 0.5

    # gas exhaust exit velocity
    return exit_velocity

def flow_rates(T_05, P_05, mach, T_01, P_01, R, gamma, exit_velocity, area_inlet, area_exit, fuel_air_ratio):
    """
    Calculate the mass flow rates through the engine

    Parameters:
        T_05 (float): Turbine outlet temperature
        P_05 (float): Turbine outlet pressure
        mach (float): Mach number
        T_01 (float): Ambient temperature (compressor inlet)
        P_01 (float): Ambient pressure (compressor inlet)
        R (float): Specific gas constant for air
        gamma (float): Specific heat ratio of the air
        exit_velocity (float): Nozzle exit velocity
        area_inlet (float): Inlet area (m²)
        area_exit (float): Nozzle exit area (m²)
        fuel_air_ratio (float): Fuel-to-air ratio

    Returns:
        tuple: Contains:
            mass_fuel_flow_rate (float): Fuel mass flow rate (kg/s)
            mass_air_flow_rate (float): Air mass flow rate (kg/s)
            mass_exit_flow_rate (float): Total exit mass flow rate (kg/s)
            inlet_velocity (float): Inlet velocity (m/s)
    """
    # calculate densities from atmospheric conditions
    rho = modules.cp.PropsSI('D', 'T', T_01, 'P', P_01, 'Air')
    rho_e = modules.cp.PropsSI('D', 'T', T_05, 'P', P_05, 'Air')

    # calculate vehicle speed in terms of m/s
    a = (gamma * R * T_01) ** 0.5
    
    # calculate fuel flow rate in kg/s
    if mach == 0:
        inlet_velocity = 55 # m/s, assuming catapult takeoff, aircraft carrier etc.
    else: 
        inlet_velocity = a * mach

    mass_air_flow_rate = rho * inlet_velocity * area_inlet
    mass_fuel_flow_rate = mass_air_flow_rate * fuel_air_ratio
    mass_exit_flow_rate = mass_fuel_flow_rate + mass_air_flow_rate

    return mass_fuel_flow_rate, mass_air_flow_rate, mass_exit_flow_rate, inlet_velocity

def thrust(mass_air_flow_rate, fuel_air_ratio, exit_velocity, inlet_velocity, P_05, P_01, area_exit, mass_fuel_flow_rate, Q_R, mach):
    """
    Calculate the thrust produced by the engine

    Parameters:
        mass_air_flow_rate (float): Air mass flow rate (kg/s)
        fuel_air_ratio (float): Fuel-to-air ratio
        exit_velocity (float): Nozzle exit velocity (m/s)
        inlet_velocity (float): Inlet velocity (m/s)
        P_05 (float): Pressure at the nozzle exit
        P_01 (float): Ambient pressure at the inlet
        area_exit (float): Nozzle exit area (m²)
        thermal_efficiency (float): Thermal efficiency
        mass_fuel_flow_rate (float): Fuel mass flow rate (kg/s)
        Q_R (float): Specific energy release of the fuel (J/kg)
        mach (float): Mach number

    Returns:
        float: Thrust produced (N)
    """
    if mach == 0:
        thermal_efficiency = ((1 + fuel_air_ratio) * (exit_velocity ** 2 / 2) - (inlet_velocity ** 2 / 2)) / (fuel_air_ratio * Q_R)
        thrust_val = 2 * thermal_efficiency * Q_R * mass_fuel_flow_rate / exit_velocity
    else: 
        # calculate thrust using conservation of momentum
        thrust_val = mass_air_flow_rate * ((1 + fuel_air_ratio) * exit_velocity - inlet_velocity) + (P_05 - P_01) * area_exit

    return thrust_val

def efficiencies(thrust, inlet_velocity, exit_velocity, fuel_air_ratio, mass_fuel_flow_rate, mass_air_flow_rate, Q_R, mach):  
    """
    Calculate the engine efficiencies and thrust specific fuel consumption (TSFC)

    Parameters:
        thrust (float): Engine thrust (N)
        inlet_velocity (float): Inlet velocity (m/s)
        exit_velocity (float): Nozzle exit velocity (m/s)
        fuel_air_ratio (float): Fuel-to-air ratio
        mass_fuel_flow_rate (float): Fuel mass flow rate (kg/s)
        mass_air_flow_rate (float): Air mass flow rate (kg/s)
        Q_R (float): Specific energy release of the fuel (J/kg)
        mach (float): Mach number

    Returns:
        tuple: Contains:
            thermal_efficiency (float): Thermal efficiency
            prop_efficiency (float): Propulsive efficiency
            overall_efficiency (float): Overall efficiency
            TSFC (float): Thrust specific fuel consumption (kg/N·hr)
    """
    # Engine efficiency properties
    thermal_efficiency = ((1 + fuel_air_ratio) * (exit_velocity ** 2 / 2) - (inlet_velocity ** 2 / 2)) / (fuel_air_ratio * Q_R)
    if mach == 0:
        prop_efficiency = 0.4
    else:
        prop_efficiency = (thrust * inlet_velocity) / (mass_air_flow_rate * ((1 + fuel_air_ratio) * 
                                                                         (exit_velocity ** 2 / 2) - (inlet_velocity ** 2 / 2)))
    overall_efficiency = thermal_efficiency * prop_efficiency

    # Thrust Specific Fuel Consumption in kg/N/hr
    TSFC = mass_fuel_flow_rate * 3600 / thrust

    return thermal_efficiency, prop_efficiency, overall_efficiency, TSFC

def aircraft_range(n_o, C_L, C_D, wet_mass, dry_mass, Q_R, thrust, mass_fuel_flow_rate, mach, T_01, P_01, wing_area):
    """
    Calculate the range of the aircraft
    
    Parameters:
        n_o (float): Overall efficiency
        C_L (float): Lift coefficient
        C_D (float): Drag coefficient
        wet_mass (float): Wet mass of the aircraft
        dry_mass (float): Dry mass of the aircraft
        Q_R (float): Specific energy release of the fuel (J/kg)
        thrust (float): Engine thrust (N)
        mass_fuel_flow_rate (float): Fuel mass flow rate (kg/s)
        mach (float): Mach number
        T_01 (float): Ambient temperature at the inlet
        P_01 (float): Ambient pressure at the inlet
        wing_area (float): Wing area (m²)
    
    Returns:
        float: Range of the aircraft
    """
    if mach == 0:
        u = 55 # m/s, assuming catapult takeoff, aircraft carrier etc.
    else: 
        u = mach * modules.CP.PropsSI('A', 'T', T_01, 'P', P_01, 'air')
    
    rho = modules.CP.PropsSI('D', 'T', T_01, 'P', P_01, 'air')

    
    lift = 0.5 * C_L * rho * u ** 2 * wing_area
    drag = 0.5 * C_D * rho * u ** 2 * wing_area

    range_1 = n_0 * (lift / drag) * modules.np.log(wet_mass / dry_mass) * (Q_R / (9.81))
    range_2 = (lift / drag) * modules.np.log(wet_mass / dry_mass) * (thrust / mass_fuel_flow_rate) * (u / 9.81)
    
    if range_1 != range_2:
        print("Warning: The two methods for calculating range are not equal.")
        print(f"Range Method 1: {range_1} km")
        print(f"Range Method 2: {range_2} km")
    else: 
        print(f"Range: {range_1} km")
        
    return range_1, range_2
    
    

# ------------------ Create Lists ------------------#
T_01 = []
P_01 = []
T_02 = []
P_02 = []
atmos_con_list = []
T_03 = []
P_03 = []
Cp = []
fuel_air_ratio = []
T_05 = []
P_05 = []
exit_velocity = []
mass_fuel_flow_rate = []
mass_air_flow_rate = []
mass_exit_flow_rate = []
inlet_velocity = []
thrust_list = []
thermal_efficiency = []
prop_efficiency = []
overall_efficiency = []
TSFC_list = []

# ------------------ Variables ------------------#
# Modelled after the -, for now, want to find some optimisation formula
# kg, dry mass, could become user defined input, doesnt effect anything for now
engine_TJ_mass = 1
max_TJ_thrust = 1  # Newtons, could become user defined input, doesnt effect anything for now
max_TJ_mach = 2  # transition speed to ramjet, to be revised
gamma = 1.4  # for air, pre combustion
pressure_ratio = 10  # for now, want this to become an user defined input eventually
max_TJ_temp = 1500  # Kelvin
Q_R = 48E6  # J/kg, C10 H22
R = 287  # J / (kg * K)
area_inlet = 0.4    # m², typical turbojet inlet area
area_exit = 0.15   # m², typical exhaust/nozzle throat area
area_ratio = area_exit / area_inlet # approximately 0.2857, want to optimise this
# Based on F14-A
C_L = 1.5  # Lift Coefficient, for now
C_D = 0.016  # Drag Coefficient, for now
wing_area = 52.5 # m², for now


# ------------------ User Defined Inputs ------------------#
mach_minimum = 0
mach_step = 0.05
mach_maximum = float(input("Enter the maximum Mach number: "))
mach_values = modules.np.arange(mach_minimum, mach_maximum + mach_step, mach_step).tolist()

altitude_minimum = 0.0183 # Height of aircraft carrier runway
altitude_maximum = float(input("Enter the maximum altitude in km: "))

if altitude_maximum == 0:
    altitude_values = [0] * len(mach_values)
else:
    altitude_step = altitude_maximum / len(mach_values)
    altitude_values = modules.np.arange(altitude_minimum, altitude_maximum, altitude_step).tolist()

# ------------------ Loop through Speeds ------------------#
for i in range(len(mach_values)):
    # inputs: gamma, mach number, altitude in km
    t_01, p_01, t_02, p_02 = compressor_inlet(gamma, mach_values[i], altitude_values[i])

    # inputs: pressure ratio (variable), compressor inlet pressure, compressor inlet temperature, atmospheric conditions, gamma
    t_03, p_03, cp = compressor_outlet(pressure_ratio, p_02, t_02, gamma, t_01, p_01)

    # inputs: maximum temperature, compressor outlet temperature, specific energy release of the fuel, compressor outlet pressure, specific heat
    f_a_ratio = fuel_to_air_ratio(max_TJ_temp, t_03, Q_R, p_03, cp)

    # inputs: maximum temperature, compressor outlet temperature, compressor inlet temperature, compressor outlet pressure
    t_05, p_05 = turbine_outlet(max_TJ_temp, t_03, t_02, p_03, gamma)

    # inputs: nozzle exit temperature, nozzle exit pressure, specific gas constant for air
    v_e = nozzle_velocity(t_05, p_05, R, p_01)

    # inputs: nozzle exit temperature, nozzle exit pressure, mach number, ambient temperature, ambient pressure, specific gas constant, gamma
    #         nozzle exit velocity, inlet area, exhaust area, fuel to air ratio
    mass_f_f_r, mass_a_f_r, mass_e_f_r, v_i = flow_rates(t_05, p_05, mach_values[i], t_01, p_01, R, gamma, v_e, area_inlet, area_exit, f_a_ratio)

    # inputs: air mass flow rate, fuel to air ratio, nozzle exit velocity, inlet velocity, nozzle exit pressure, ambient pressure, nozzle area
    thrust_val = thrust(mass_a_f_r, f_a_ratio, v_e, v_i, p_05, p_01, area_exit, mass_f_f_r, Q_R, mach_values[i])

    # inputs: thrust, inlet velocity, exhaust velocity, fuel to air ratio, fuel mass flow rate, air mass flow rate, specific energy release of the fuel
    n_th, n_p, n_o, TSFC_val = efficiencies(thrust_val, v_i, v_e, f_a_ratio, mass_f_f_r, mass_a_f_r, Q_R, mach_values[i])
    
    # 
    range_1, range_2 = aircraft_range(n_o, C_L, C_D, engine_TJ_mass, engine_TJ_massm, mass_f_f_r, Q_R, thrust_val, mass_f_f_r, mach    

    # Append values to lists
    T_01.append(t_01)
    P_01.append(p_01)
    T_02.append(t_02)
    P_02.append(p_02)
    T_03.append(t_03)
    P_03.append(p_03)
    Cp.append(cp)
    fuel_air_ratio.append(f_a_ratio)
    T_05.append(t_05)
    P_05.append(p_05)
    exit_velocity.append(v_e)
    mass_fuel_flow_rate.append(mass_f_f_r)
    mass_air_flow_rate.append(mass_a_f_r)
    mass_exit_flow_rate.append(mass_e_f_r)
    inlet_velocity.append(v_i)
    thrust_list.append(thrust_val)
    thermal_efficiency.append(n_th)
    prop_efficiency.append(n_p)
    overall_efficiency.append(n_o)
    TSFC_list.append(TSFC_val)


# ------------------ Transform data for Plotting ------------------#
# Temperatures and Pressures for different stages of T-S, P-V diagrams
temperatures = [modules.stat.mean(T_01), modules.stat.mean(
    T_02), modules.stat.mean(T_03), max_TJ_temp, modules.stat.mean(T_05)]
pressure = [modules.stat.mean(P_01), modules.stat.mean(P_02), modules.stat.mean(
    P_03), modules.stat.mean(P_03), modules.stat.mean(P_05)]

# Calculate specific volume and entropy
specific_volumes = [R * T / P for T, P in zip(temperatures, pressure)]
entropies = [modules.cp.PropsSI('S', 'T', T, 'P', P, 'Air')
             for T, P in zip(temperatures, pressure)]

pressure_kpa = [P / 1000 for P in pressure]
thrust_array = modules.np.array(thrust_list)

# ------------------ Plotting ------------------#
modules.plt.figure()
modules.plt.plot(entropies, temperatures)
modules.plt.plot([entropies[0], entropies[-1]], [temperatures[0], temperatures[-1]],
                 ':', color='gray', linewidth=1, label='Start/End Connection')
modules.plt.xlabel('Entropy (J/kg·K)')
modules.plt.ylabel('Temperature (K)')
modules.plt.title('T-S Diagram')

modules.plt.figure()
modules.plt.plot(specific_volumes, pressure_kpa)
modules.plt.plot([specific_volumes[0], specific_volumes[-1]], [pressure_kpa[0], pressure_kpa[-1]],  # Connect start and end
                 ':', color='gray', linewidth=1, label='Start/End Connection')
modules.plt.xlabel('Specific Volume (m³/kg)')
modules.plt.ylabel('Pressure (kPa)')
modules.plt.title('P-V Diagram')

modules.plt.figure()
modules.plt.plot(mach_values, thrust_array)
modules.plt.xlabel('Mach Number')
modules.plt.ylabel('Thrust (kN)')
modules.plt.title('Thrust vs Mach Number')

modules.plt.figure()
modules.plt.plot(altitude_values, thrust_array)
modules.plt.xlabel('Altitude (km)')
modules.plt.ylabel('Thrust (kN)')
modules.plt.title('Thrust vs Altitude')

modules.plt.show()

# ------------------ Print important values ------------------#
print(f"The turbojet engine produces {max(thrust_list)/1e3:.2f}kN of Peak Thrust at Mach {mach_values[thrust_list.index(max(thrust_list))]:.3f}")
print(f"Thermal Efficiency: {modules.stat.mean(thermal_efficiency):.3f} \nPropulsive Efficiency: {modules.stat.mean(prop_efficiency):.3f} \
      \nOverall Efficiency: {modules.stat.mean(overall_efficiency):.3f}")
print(f"Peak Overall Efficiency: {max(overall_efficiency):.3f} at Mach {mach_values[overall_efficiency.index(max(overall_efficiency))]:.3f}")
print(f"The engine consumes, on average, {modules.stat.mean(mass_fuel_flow_rate):.4f} kg of fuel per second")
print(f"Thrust Specific Fuel Consumption: {modules.stat.mean(TSFC_list):.4f} kg/N·hr")
print(f"At 0 Mach, the engine produces {thrust_list[0]/1e3:.2f}kN of Peak Thrust")

