using Pkg
using Ipopt, JuMP
using Dates
using CSV
using DataFrames

println("--- Start Program ---")

# GENERAL PARAMETERS
tfinal = 10;
dt = 1;                     # (h) timestep
years=15;
efficiency_bat = 0.97;      # check with professor
eta_charge = 0.95;          # check with professor
eta_discharge = 0.95;       # check with professor
SOC_bat_MAX = 1;            # (-) - Maximum SOC for batteries
SOC_bat_MIN = 0.35;         # (-) - Minimum SOC for batteries

# COST PARAMETERS
diesel_capex = 823e3;                       # Capex diesel
diesel_opex = 150 *years;                   # Opex diesel
solar_capex = 1065e3;                       # Capex Solar 
solar_opex = 14800*years;                   # Opex Solar
bat_capex = 600*1000;                       # Capex Battery
bat_opex = 7240*years;                      # Opex Battery

# DEMAND PARAMETERS
demand_t = CSV.read("data/Cornell_heating_load_2017.csv", DataFrame);

# GENERATION PARAMETERS
solar_available = CSV.read("data/PV_Gen_CatalunyA_1KW.csv", DataFrame);

## INITIAL CONDITIONS
SOC_ini = 0.6;                     # (p.u) - Initial state of charge of the battery

# Model

#model
m = Model(Ipopt.Optimizer)

## VARIABLES
@variable(m, diesel_generation_t[1:tfinal] >= 0)  # (MW) diesel gen for each time step t
@variable(m, diesel_capacity >= 0)  # (MWh) diesel capacity
@variable(m, solar_generation_t[1:tfinal] >= 0)  # (MW) solar gen for each time step t
@variable(m, solar_capacity >= 0)  # (MWh) solar capacity
@variable(m, battery_energy_capacity >= 0)  # (MWh) battery capacity
@variable(m, battery_power_capacity >= 0)  # (MW) battery capacity
@variable(m, charge_battery_t[1:tfinal] >= 0)  # (MW) - Charge power for the battery
@variable(m, discharge_battery_t[1:tfinal] >= 0)  # (MW) - Discharge power for the battery
@variable(m, SOC_battery[1:tfinal] >= 0)  # (p.u) - State of charge of the battery 

# OBJECTIVE FUNCTION

@objective(m, Min, diesel_capex*diesel_capacity + sum(diesel_opex*diesel_generation_t[1:tfinal]) + sum(solar_opex*solar_generation_t[1:tfinal]) + solar_capex + bat_opex*battery_energy_capacity + bat_capex);

# CONSTRAINT 1: DIESEL GENERATION FOR ANY HOUR MUST BE LESS THAN MAX CAPACITY
for ti = 1:tfinal
    @NLconstraint(m, diesel_generation_t[ti] <= diesel_capacity);
end

# CONSTRAINT 2: SOLAR GENERATION FOR ANY HOUR MUST BE LESS THAN MAX CAPACITY
for ti = 1:tfinal
    @NLconstraint(m, solar_generation_t[ti] <= solar_available[ti,3]);
end

# CONTRAINTS 3: BATTERY CHARGE FOR ANY HOUR MUST BE LESS THAN MAX
for ti = 1:tfinal
    @NLconstraint(m, charge_battery_t[ti] <= battery_power_capacity);
end

# COSTRAINT 4: BATTERY DISCHARGE FOR ANY HOUR MUST BE LESS THAN MAX
for ti = 1:tfinal
    @NLconstraint(m, discharge_battery_t[ti] <= battery_power_capacity);
end

# CONSTRAINT 5: DISCHARGE CAPACITY IS HALF THE BATTERY POWER CAPACITY
@NLconstraint(m, battery_power_capacity == 0.5*battery_energy_capacity);

# CONSTRAINTS 6: STATE OF CHARGE TRACKING
@NLconstraint(m, SOC_battery[1] == SOC_ini);

for ti = 2:tfinal
    @NLconstraint(m, SOC_battery[ti] == SOC_battery[ti-1] + (((eta_charge*charge_battery_t[ti])-(discharge_battery_t[ti]/eta_discharge))*dt)/battery_energy_capacity);
end

# CONSTRAINT 7: DEMAND BALANCING
for ti = 1:tfinal
    @NLconstraint(m, diesel_generation_t[ti] + solar_capacity*solar_generation_t[ti] - charge_battery_t[ti] + discharge_battery_t[ti] == demand_t[ti,2]);
end

# CONSTRAINT 8a: SOC LIMITS (MAXIMUM)
for ti = 1:tfinal
    @NLconstraint(m, SOC_bat_MAX >= SOC_battery[ti]);
end

# CONSTRAINT 8b: SOC LIMITS (MINIMUM)
for ti = 1:tfinal
    @NLconstraint(m, SOC_battery[ti] >= SOC_bat_MIN);
end

## EXECUTION OF OPTIMIZATION PROBLEM
optimize!(m)


#### EXPORT DATA TO CSV #######

#Store values hourly

#Diesel
diesel_cap_opt = JuMP.value.(diesel_capacity)
diesel_cap_opt_list = zeros(tfinal)
diesel_cap_opt_list[1] = diesel_cap_opt
diesel_gen_opt = JuMP.value.(diesel_generation_t)

#diesel_opt = DataFrame(Diesel_Capacity_MW = diesel_cap_opt_list, Diesel_generation_in_hour=diesel_gen_opt)

#solar
solar_cap_opt = JuMP.value.(solar_capacity)
solar_cap_opt_list = zeros(tfinal)
solar_cap_opt_list[1] = solar_cap_opt
solar_gen_opt = JuMP.value.(solar_generation_t)

#solar_opt = DataFrame(Solar_Capacity_MW = solar_cap_opt_list, Solar_generation_in_hour=solar_gen_opt)

#battery

batt_Ecap_opt = JuMP.value.(battery_energy_capacity)
batt_Ecap_opt_list = zeros(tfinal)
batt_Ecap_opt_list[1] = batt_Ecap_opt
batt_Pcap_opt = JuMP.value.(battery_power_capacity)
batt_Pcap_opt_list = zeros(tfinal)
batt_Ecap_opt_list[1] = batt_Pcap_opt
batt_charge_opt = JuMP.value.(charge_battery_t)
batt_discharge_opt = JuMP.value.(discharge_battery_t)
batt_SOC_opt = JuMP.value.(SOC_battery)

#batt_opt = DataFrame(Battery_Energy_Cap_MWh = batt_Ecap_opt_list,  Battery_Power_Cap_MWh = batt_Pcap_opt_list, Battery_Charge_Cap_MW =batt_charge_opt, Battery_Disharge_Cap_MW =batt_discharge_opt, Battery_SOC =  batt_SOC_opt)

overall_opt = DataFrame(hour= 1:tfinal,Diesel_Capacity_MW = diesel_cap_opt_list, Diesel_generation_in_hour=diesel_gen_opt,Solar_Capacity_MW = solar_cap_opt_list, Solar_generation_in_hour=solar_gen_opt,Battery_Energy_Cap_MWh = batt_Ecap_opt_list,  Battery_Power_Cap_MWh = batt_Pcap_opt_list, Battery_Charge_Cap_MW =batt_charge_opt, Battery_Disharge_Cap_MW =batt_discharge_opt, Battery_SOC =  batt_SOC_opt)

CSV.write("results/Optimal_Values.csv", overall_opt)


##### CHECK DATA RESULTS ON CONSOL #####

println("Diesel Cap: ")
println(JuMP.value.(diesel_capacity));
println("Battery Energy Cap: ")
println(JuMP.value.(battery_energy_capacity));