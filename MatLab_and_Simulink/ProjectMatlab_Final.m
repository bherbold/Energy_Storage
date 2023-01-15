%% INPUT DATA ----------------------------------------------------------- %
clear all
clc

% Read CSV Data
PV_Capacity = 315.85; %KW
pv_gen = readtable('Data/PV_Germany_2019');

criticalh = 7048;%most critical hour

starth = 4057; %first hour % % criticalh
endh= 4080; %second hour criticalh +1%4080;%criticalh + 1;
secinh = 3600; %in how many values should an hour be divided
totalstepnum = (1+endh-starth)*secinh
    
% Base system values
wbase = 2*pi*50; % [rad/s] - Base electrical frequency
Psysbase = 320; % [kW] - Base system power
% Load
load_0 = 0.5; % [p.u.] - Initial load w.r.t. base system power
load_f = 0.7; % [p.u.] - Final load w.r.t. base system power
% Renewables
penetration = 0.8; % [p.u.] - Renewables penetration w.r.t. base system power
% Synchronous Generators
Pgenbase = Psysbase; % [MW] - Rated power
TG = 0.2; % [-] - Turbine parameter. Source: (Kundur, pg 599)
TCH = 0.3; % [-] - Turbine parameter. Source: (Kundur, pg 599)
R = 0.03; % [-] - Speed governor droop characteristic
H = 5; % [s] - Inertia constant
k = 0.5; 

%pv_gen = pv_gen{:,3};4057 4080
pv_array = [];

for c = starth:endh
    help = (pv_gen{c,2} * PV_Capacity)/Psysbase; %in kW
    for i= 1:secinh
        pv_array = [pv_array, help];
    end
end
pv_array = pv_array';
pv_data = timeseries(pv_array);
%ttt.TimeInfo.Units = 'hours';
pv_data;

load_profile = readtable('Data/Hourly_Load_2018_Tool_Manufacturer.xlsx');

load_array = [];

for c = starth:endh
    help = load_profile{c,2}/Psysbase; % in kW
    for i= 1:secinh
        load_array = [load_array, help];
    end
end

load_array = load_array';
timeseries_load = timeseries(load_array);

%Diesel Generator Capacity
diesel_cap = 91;%42.43;
diesel_cap_pu = diesel_cap/Psysbase;

%find highest load difference




%% BATTERY-BASED ENERGY STORAGE SOLUTION (LITHIUM-ION)
% RATED POWER, ENERGY AND VOLTAGE

bat_power_cap = 106.98; % power capacity in kw
bat_energy_cap = 755.03; % power capacity in kw

storage_p = bat_power_cap/Psysbase; % [p.u.] - Storage rated power output w.r.t. base system power
storage_e = bat_energy_cap/Psysbase; % [p.u.] - Storage rated energy capacity (2 is the number of hours
%storage_capacity =  0.7 * Psysbase;%MW
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
EFF = 0.93; % (-) - Efficiency
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
Rf = 0.2;
pdboost = 0.7;
k_bat = 0.5;
load_0_bat = (load_0 - penetration)*Psysbase/Pgenbase; 
load_0_synch = 0;%(load_0 - penetration )*Psysbase/Pgenbase; ;%(load_0 - penetration )*Psysbase/Pgenbase; % [p.u.] - Initial setpoint for synchronous gen.

R_PV = 0.1;
k_PV = 0.5;

%% INITIAL CONDITIONS FOR SIMULATION
SOC_0 = 0.5;
OCV_0 = interp1(SOC_LUT,OCV_LUT,SOC_0,'linear','extrap');

%% Initial Values
best_setting = [];
Rf_array = [];

out = sim('ProjectSimulink.slx');
magnitude = max(out.w.Data);
Rf_array = [Rf_array; Rf magnitude];
[best_Rf_M, best_Rf_index] = min(Rf_array(:,2));
Rf = Rf_array(best_Rf_index,1);
best_setting = [best_setting; Rf_array(best_Rf_index,1),best_Rf_M];


%% Sensitivity Analysis
counter1 = 0;
best_setting = [];

Rf_array = [];
for Rf = 0.05:0.05:0.25
out = sim('ProjectSimulink.slx');
magnitude = max(out.w.Data);
Rf_array = [Rf_array; Rf magnitude];

end

[best_Rf_M, best_Rf_index] = min(Rf_array(:,2));
Rf = Rf_array(best_Rf_index,1);
best_setting = [best_setting; Rf_array(best_Rf_index,1),best_Rf_M];

kbat_array = [];
for k_bat = 0.5:0.1:0.9
out = sim('ProjectSimulink.slx');
magnitude = max(out.w.Data);
kbat_array = [kbat_array; k_bat magnitude];

end
[best_kbat_M, best_kbat_index] = min(kbat_array(:,2));
k_bat = kbat_array(best_kbat_index,1);
best_setting = [best_setting; kbat_array(best_kbat_index,1),best_kbat_M];

pdboost_array = [];
for pdboost = 0.6:0.1:1
out = sim('ProjectSimulink.slx');
magnitude = max(out.w.Data);
pdboost_array = [pdboost_array; pdboost magnitude];
conuter1 = counter1 +1
end
[best_pdboost_M, best_pdboost_index] = min(pdboost_array(:,2));
pdboost = pdboost_array(best_pdboost_index,1);
best_setting = [best_setting; pdboost_array(best_pdboost_index,1),best_pdboost_M];

R_array = [];
for R = 0.020:0.005:0.035
out = sim('ProjectSimulink.slx');
magnitude = max(out.w.Data);
R_array = [R_array; R magnitude];
conuter1 = counter1 +1
end
[best_R_M, best_R_index] = min(R_array(:,2));
R = R_array(best_R_index,1);
best_setting = [best_setting; R_array(best_R_index,1),best_R_M];
%%
k_array = [];
for k = 0.5:0.1:0.6
out = sim('ProjectSimulink.slx');
magnitude = max(out.w.Data);
k_array = [k_array; k magnitude];
conuter1 = counter1 +1
end
[best_k_M, best_k_index] = min(k_array(:,2));
k = k_array(best_k_index,1);
best_setting = [best_setting; k_array(best_k_index,1),best_k_M];


%% Sensitivity PV 

% Only makes sense when SOC over 80%

SOC_0 = 0.95;
OCV_0 = interp1(SOC_LUT,OCV_LUT,SOC_0,'linear','extrap');

R_PV_array = [];
for R_PV = 0.1:0.1:0.5
out = sim('ProjectSimulink.slx');
magnitude = min(out.w.Data);
R_PV_array = [R_PV_array; R_PV magnitude];
conuter1 = counter1 +1
end
[best_R_PV_M, best_R_PV_index] = max(R_PV_array(:,2));
R_PV = R_PV_array(best_R_PV_index,1);
best_setting = [best_setting; R_PV_array(best_R_PV_index,1),best_R_PV_M];
%%
k_PV_array = [];
for k_PV = 0.5:0.1:0.9
out = sim('ProjectSimulink.slx');
magnitude = min(out.w.Data);
k_PV_array = [k_PV_array; k_PV magnitude];
conuter1 = counter1 +1
end
[best_k_PV_M, best_k_PV_index] = max(k_PV_array(:,2));
k_PV_array = [k_PV_array; k_PV magnitude];
best_setting = [best_setting; k_PV_array(best_k_PV_index,1),best_k_PV_M];

%print(best_setting);
%% Final Model with Generator
SOC_0 = 0.7;
OCV_0 = interp1(SOC_LUT,OCV_LUT,SOC_0,'linear','extrap');

Rf = 0.05;
k_bat = 0.5;
pdboost = 1;
R = 0.02;
k = 0.5;
R_PV = 0.1;
k_PV = 0.5;

%out = sim('ProjectSimulink.slx');



%% Execute model
figure(2)
title('Frequency (Hz)')
hold on
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ylim([49.9 50.2])
colors = {'r','b','g','c'};
grid on
box on

plot(out.w.Time(100000:180000),out.w.Data(100000:180000),'linewidth',2)
%%

figure(1)
title('Generation')
hold on
xlabel('Time (s)')
ylabel('KW')
ylim([-40 80])
colors = {'r','b','g','c'};
grid on
box on

plot(out.Pmech.Time(100000:180000),out.Pmech.Data(100000:180000),'linewidth',2)
hold on
plot(out.P_bat.Time(100000:180000),out.P_bat.Data(100000:180000),'linewidth',2)
hold on
plot(out.PV_gen.Time(100000:180000),out.PV_gen.Data(100000:180000),'linewidth',2)
legend('Diesel','Battery','PV','Location','Best')

%%

figure(1)
title('SOC')
hold on
xlabel('Time (s)')
ylabel('%')
ylim([0 1])
colors = {'r','b','g','c'};
grid on
box on

plot(out.SOC.Time,out.SOC.Data,'linewidth',2)


%%

i = 1;


for Rf = 0.1:0.05:0.3
out = sim('ProjectSimulink.slx');
plot(out.w.Time,out.w.Data,'linewidth',2)
i = i+1;
end
legend('Rf = 0.1','Rf = 0.15','Rf = 0.2','Rf = 0.25','Rf = 0.3','Location','Best')
% -----------------------------------------------------------------------

%plot(out.w.Time,out.w.Data,'linewidth',2);
%%
R = 0.03;
out = sim('ProjectSimulink.slx');

%% Execute model
figure(2)
title('Changing R Diesel Value - Frequency')
hold on
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ylim([49 50.2])
colors = {'r','b','g','c'};
grid on
box on



i = 1;


for k = 0.3:0.05:0.45
out = sim('ProjectSimulink.slx');
plot(out.w.Time,out.w.Data,'linewidth',2)
i = i+1;
end
legend('R = 0.3','R = 0.35','R = 0.4','R = 0.45','Location','Best')
% -----------------------------------------------------------------------


%% Graphics

% figure(2)
% title('Changing Rf Battery Value - Generation')
% hold on
% xlabel('Time (s)')
% ylabel('KW')
% ylim([49 50.2])
% colors = {'r','b','g','c'};
% grid on
% box on
% 
% plot(out.Pmech.Time,(out.Pmech.Data*Psysbase),'linewidth',2)
% hold on
% plot(out.P_bat.Time,(out.P_bat.Data*Psysbase),'linewidth',2);
% hold on
% plot(out.PV_gen.Time,(out.PV_gen.Data*Psysbase),'linewidth',2)
% hold on
% plot(out.demand.Time,(out.demand.Data*Psysbase),'linewidth',2);

%SOC_chart = plot(out.SOC.Time,out.SOC.Data,'linewidth',2)
%freq_chart = plot(out.w.Time,out.w.Data,'linewidth',2);



%exportgraphics(SOC_chart,'Graphs/SOC_Chart.png','Resolution',300);
%hold on
%plot(out.results.time,out.Battery_results.signals(2).values,'linewidth',2)
