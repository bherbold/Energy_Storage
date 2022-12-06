%% INPUT DATA ----------------------------------------------------------- %
clear all
clc

% Read CSV Data
PV_Capacity = 100; %KW
pv_gen = readtable('Data/PV_Gen_CatalunyA_1KW.csv');
%pv_gen = pv_gen{:,3};
pv_array = [];

for c = 1:8760
    help = pv_gen{c,3} * PV_Capacity/1000; %in MW
    %for i= 1:3600
        pv_array = [pv_array, help];
    %end
end
pv_array = pv_array';
ttt = timeseries(pv_array);
%ttt.TimeInfo.Units = 'hours';
ttt
    
% Base system values
wbase = 2*pi*50; % [rad/s] - Base electrical frequency
Psysbase = 100; % [MW] - Base system power
% Load
load_0 = 0.7; % [p.u.] - Initial load w.r.t. base system power
load_f = 0.5; % [p.u.] - Final load w.r.t. base system power
% Renewables
penetration = 0.3; % [p.u.] - Renewables penetration w.r.t. base system power
% Synchronous Generators
Pgenbase = Psysbase; % [MW] - Rated power
TG = 0.2; % [-] - Turbine parameter. Source: (Kundur, pg 599)
TCH = 0.3; % [-] - Turbine parameter. Source: (Kundur, pg 599)
R = 0.03; % [-] - Speed governor droop characteristic
H = 5; % [s] - Inertia constant
k = 0.1; 
load_0_synch = (load_0 - penetration)*Psysbase/Pgenbase; % [p.u.] - Initial setpoint for synchronous gen.

%% BATTERY-BASED ENERGY STORAGE SOLUTION (LITHIUM-ION)
% RATED POWER, ENERGY AND VOLTAGE
storage_p = 1; % [p.u.] - Storage rated power output w.r.t. base system power
storage_e = storage_p * 2; % [p.u.] - Storage rated energy capacity (2 is the number of hours
% we want to get the rated power storage_p
storage_u = 600; % [V] - Rated voltage of the battery pack
% CELL PARAMETERS
R_cell = 0.01; % (Ohm) - Discharge resistance
Cn_cell = 7200; % (As) - Cell capacity
OCV_LUT_cell = [3.2472 3.4658 3.5546 3.5987 3.6254 3.6645 3.7531 3.8397 3.9400 4.0502 4.1763]';
% PACK PARAMETERS
Ns = ceil( storage_u / mean(OCV_LUT_cell) ); % (-) - Number of cells in series
Np = ceil( storage_e * Psysbase * 1e6 / ( Ns * mean(OCV_LUT_cell) * Cn_cell / 3600 ) ); % (-) - Number of cells in parallel
Rb = R_cell*Ns/Np; % (Ohm) - Pack resistance
Cn = Cn_cell*Np; % (As) - Pack capacitance
OCV_LUT = Ns*OCV_LUT_cell';% (V) - Pack OCV
SOC_LUT = [0.0066 0.1004 0.2004 0.3003 0.4002 0.5002 0.6001 0.7000 0.8000 0.8999 0.9998]';
Ub_rated = mean(OCV_LUT); % (V) - Pack rated OCV
Ub_exp = max(OCV_LUT); % (V) - Pack maximum OCV
Ub_min = min(OCV_LUT); % (V) - Pack minimum OCV
EFF = 0.95; % (-) - Efficiency
% DC-DC CONVERTER PARAMETERS
udc = storage_u * 1.5; % (V) - DC-link voltage
L = 0.005; % (H) - Inductance between converter and battery
Irated = storage_p * Psysbase * 1e6 / Ub_min; % (A) - Maximum admissible current (battery side)
% DC-DC CONVERTER - CURRENT CONTROLLER
% TO BE SELECTED BY THE DESIGNER
wi = 10; % (rad/s) - Natural frequency
xii = 0.7; % (-) - Damping coefficient
% CALCULATIONS
Kpi = 2*xii*wi*L/udc; % (-) - Proportional control gain
Kii = wi^2*L/udc; % (-) - Integral control gain
% DC-DC CONVERTER - SOC CONTROLLER
Ksoc = 100*(Irated*0.05);
% END USER CONTROL ALGORITHM - FREQUENCY CONTROL
Rf = 0.03;

%% INITIAL CONDITIONS FOR SIMULATION
SOC_0 = 0.7;
OCV_0 = interp1(SOC_LUT,OCV_LUT,SOC_0,'linear','extrap');

%% Execute model

out=sim('ProjectSimulink.slx')

%plot(out.results.time,out.Battery_results.signals(2).values,'linewidth',2)
%hold on
%plot(out.results.time,out.Battery_results.signals(2).values,'linewidth',2)
