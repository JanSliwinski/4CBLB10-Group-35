warning off
%% Initialization
if ~exist("T","var")
    load('T.mat','T')
end
clearvars -except T; close all; clc
addpath("Functions", "Nasa");  % Add necessary paths
% figure('Visible', 'on');
%% Units and Constants
mm = 1e-3;
dm = 0.1;
bara = 1e5;
MJ = 1e6;
kWhr = 1000 * 3600;
volperc = 0.01;  % Emissions are in volume percentages
ppm = 1e-6;      % Some are in ppm (also a volume fraction)
g = 1e-3;
s = 1;
RPM = 1500;     % constant RPM of experiments [rotation/min]
T_int = 295.15; % assumed temperature at the intake [K]
M_diesel = 167; % Molar mass (g/mol)
M_HVO = 226; % Molar mass (g/mol)
Density_HVO = 0.78;
Density_Diesel = 0.85;
LHV_diesel = 43e3;  %Lower heating value given in the project guide for Diesel B7 J/g
LHV_HVO = 44.29e3;  %Lower heating value given in the project guide for HVO J/g
LHV_GTL = 44e3;     %Lower heating value given in the project guide for GTL B7 J/g
x_HVO = 16;         %carbon atoms in HVO
x_diesel = 12;      %carbon atoms in diesel 
x_GTL = 14;         %carbon atoms in GTL

%% Load and Reshape data (also exhaust data)
ID = 'L50I14C0FuelHVO'; %DEFINE ID OF THE EXPERIMENT DATA YOU WANT TO LOAD IN!
%run fucntion to load in all relevant data
[dataIn, ExhaustData, Ca, p_filt, p_avg, S_current, mfr_fuel, CO_percent_load, HC_ppm_load, NOx_ppm_load, CO2_percent_load, O2_percent_load, lambda_load] = loadingfromT(T, ID, bara);
IDsforKPI = table({'L50I14C0FuelHVO', 'L50I15C0FuelHVO'}');

%% Define Fuel used and applicable LHV - CHANGE EVERY LINE IN THIS SECTION IF RUN WITH A DIFFERENT FUEL!!!
fuel_used = 'Diesel';
perc_blend = 0; %fraction of the blended in fuel (HVO or GTL)
x_blend = x_diesel;  %carbon atoms in the given fuel, can be: x_diesel, x_HVO or x_GTL
LHV_blend = LHV_diesel;

%% Calculate LHV and x for the fuel used
perc_diesel = 1-perc_blend;
LHV = LHV_blend * perc_blend + LHV_diesel * perc_diesel;
x = x_blend * perc_blend + x_diesel * perc_diesel;

%% Load NASA Data (if needed)
global Runiv
Runiv = 8.314; 

[SpS, El] = myload('Nasa/NasaThermalDatabase.mat', {'N2', 'O2', 'CO2', 'H2O', 'Diesel'});

%% Volume
% Engine Geometry Parameters
Cyl.Bore = 104 * mm;               % Cylinder bore
Cyl.Stroke = 85 * mm;              % Cylinder stroke
Cyl.CompressionRatio = 21.5;       % Compression ratio
Cyl.ConRod = 136.5 * mm;           % Connecting rod length
Cyl.TDCangle = 180;                % Top Dead Center angle
% Calculate cylinder volume using CylinderVolume function
volume = CylinderVolume(Ca,Cyl);
disp('Cylider volume calculated / cycle');

%% Stoichiometric calculations
fuel_name = 'Diesel';
[stoich_coeffs, reaction_eq, AFR_stoich] = StoichiometricCombustion(fuel_name, SpS, El);
stoich_coeffs.fuel

%% Detect Start and End of Injection from Sensor Current

% Average sensor current across all cycles
S_current_avg = mean(S_current, 2);

% Define a baseline current
baseline_current = min(S_current_avg);

% Define a threshold above the baseline to detect injection
threshold = baseline_current + 0.2;

% Find indices where the sensor current exceeds the threshold
injection_indices = find(S_current_avg > threshold);

% Determine the start and end indices of injection
injection_start_idx = injection_indices(1);
injection_end_idx = injection_indices(end);

% Find the corresponding crank angles
injection_start_ca = Ca(injection_start_idx, 1);
injection_end_ca = Ca(injection_end_idx, 1);

% Display the results
fprintf('Injection starts at %.2f° CA\n', injection_start_ca);
fprintf('Injection ends at %.2f° CA\n', injection_end_ca);

%% Plot Average Sensor Current with Injection Markers
figure;

% Plot the average sensor current
plot(Ca(:, 1), S_current_avg, 'LineWidth', 1.5);
xlabel('Crank Angle (°)');
ylabel('Sensor Current (A)');
title('Average Sensor Current vs. Crank Angle');
grid on;
xlim([-25, 25]);

% Add markers for injection start and end
hold on;
plot(injection_start_ca, S_current_avg(injection_start_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
plot(injection_end_ca, S_current_avg(injection_end_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);

% Annotate the start and end points
text(injection_start_ca, S_current_avg(injection_start_idx), sprintf('Start: %.2f°', injection_start_ca), ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
text(injection_end_ca, S_current_avg(injection_end_idx), sprintf('End: %.2f°', injection_end_ca), ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);

hold off;

true_mfr_fuel = CalculateMassFlowFuel(mfr_fuel, S_current, Ca, RPM, threshold);

%% Plot Pressure vs. Crank Angle for All Cycles
figure;
plot(Ca, p_filt / bara, 'LineWidth', 1);
xlabel('Crank Angle (°)');
ylabel('Pressure (bar)');
xlim([-360, 360]);
ylim([0, 50]);
title('Pressure vs. Crank Angle for All Cycles');
grid on;

% Highlight a specific cycle
iselect = 10;
hold on;
plot(Ca(iselect), p_filt(iselect) / bara, 'r', 'LineWidth', 2);

% Plot valve events
% YLIM = ylim;
% line([ValveEvents.CaIVC, ValveEvents.CaIVC], YLIM, 'Color', 'b', 'LineWidth', 1);
% line([ValveEvents.CaEVO, ValveEvents.CaEVO], YLIM, 'Color', 'r', 'LineWidth', 1);

set(gca, 'XTick', -360:60:360);
grid on;

% %% Calculate Average Volume and Pressure
% V_avg = mean(volume, 2);         % Average volume across all cycles for every CA
% p_avg = mean(p_filt, 2);             % Average pressure across all cycles for every CA
% p_filtered_avg = mean(p, 2);






%% Calculate Work
W = trapz(volume, p_avg*1e5); % Calculate the area under the averaged p-V curve
disp(['Calculated work: ', num2str(W), ' J']);

%% Calculate mass flow of fuel
% Constants
M_C = 12; % Molar mass of Carbon (g/mol)
M_H = 1; % Molar mass of Hydrogen (g/mol)
M_CO2 = 44; % Molar mass of CO2 (g/mol)
mass_CH_ratio = 5.49; % Typical C/H mass ratio for diesel

% Inputs: Known CO2 mass flow rate
CO2_mass_flow_rate = 0.5; % IDK what this value is but Barrt Sommers says we should have it (we have the percentage of it but IDK how to get mfr from that)

% Convert mass CH ratio to molar CH ratio
molar_CH_ratio = mass_CH_ratio * (M_H / M_C);

% Determine y for the given fuel (x is given above) 
y = molar_CH_ratio * x;

% Calculate molar mass of diesel (CxHy)
fuel_molar_mass = x * M_C + y * M_H; % in g/mol

% Convert CO2 mass flow rate to molar flow rate
CO2_molar_flow_rate = CO2_mass_flow_rate / M_CO2; % in mol/s

% Determine the molar flow rate of diesel
fuel_molar_flow_rate = CO2_molar_flow_rate / x;

% Convert diesel molar flow rate to mass flow rate
fuel_mass_flow_rate = fuel_molar_flow_rate * fuel_molar_mass; % in g/s

% Result
fprintf('Fuel mass flow rate for diesel: %.6f g/s\n', fuel_mass_flow_rate);


%% Stoichiometric calculations for diesel
fuel_name = 'Diesel';
[stoich_coeffs, reaction_eq, AFR_stoich] = StoichiometricCombustion(fuel_name, SpS, El);
stoich_coeffs.fuel

%% Stoichiometric Calculations for HVO
m_fuel_per100cycles = 0.128; %(g)

%for HVO100
moles_HVO = m_fuel_per100cycles / M_HVO;
moles_O2 = moles_HVO * (3*x_HVO*0.5) + 1/2; % Moles of O2 required
moles_N2 = moles_O2 * 79/21;   % Moles of N2 in reactants
moles_CO2 = moles_HVO * 16; % Moles of CO2 produced
moles_H2O = moles_HVO * 17; % Moles of H2O produced
Stoich_HVO100 = [moles_HVO, moles_O2, moles_N2, moles_CO2, moles_H2O, moles_N2];

%for HVO50
 a = m_fuel_per100cycles / (((Density_HVO / Density_Diesel) + 1) * M_diesel);
 b = m_fuel_per100cycles / (((Density_Diesel / Density_HVO) + 1) * M_HVO);
 e = a * 12 + b * 16;
 g = a * 23 + b * 34;
 q = 9.25 * a + 24.5 * b;
 k = q * (79 / 21);
Stoich_HVO50 = [a, b, q, k, e, g/2, k]; %Moles of each Component per 100 cycles


%% aROHR

%p_filt = sgolayfilt(p,2,101);

%p_filtRough = sgolayfilt(p_filt,2,101);

if exist('O2_percent_load','var')
    gamma = CalculateGamma(SpS,volume,p_avg,O2_percent_load,CO2_percent_load,true_mfr_fuel,AFR_stoich,RPM);
else
    gamma = 1.32; % constant if no exhuast data exist
end
aROHR = get_aROHR(p_avg,volume,gamma);
% Isolate peak
idxStart = 355 / 0.2; idxEnd = (50 + 360) / 0.2;
aROHR(1:idxStart) = 0; aROHR(idxEnd:end) = 0;
aHR = get_aHR(aROHR);
disp("aROHR and aHR computed successfully!")


 figure;
 subplot(1,2,1)
 plot(Ca,aROHR);xlabel("Crank Angle [deg]");ylabel("Apparent Rate of Heat Realease [J/deg]");
 xlim([-10,30]);
 legend("aROHR","Location","southeast");title(["Apparent Rate of Heat", "Release for " + ID]);
 subplot(1,2,2)
 plot(Ca,aHR);xlabel("Crank Angle [deg]");ylabel("Apparent Heat Realease [J]");
 xlim([-10,30]);
 legend("aHR","Location","southeast");title(["Apparent Heat Release", "for " + ID]);

figure;
subplot(1,2,1)
plot(Ca,aROHR);xlabel("Crank Angle [deg]");ylabel("Apparent Rate of Heat Realease [J/deg]");
xlim([-5,50]);
legend("aROHR","Location","southeast");title(["Apparent Rate of Heat", "Release for " + ID]);
subplot(1,2,2)
plot(Ca,aHR);xlabel("Crank Angle [deg]");ylabel("Apparent Heat Realease [J]");
xlim([-5,50]);
legend("aHR","Location","southeast");title(["Apparent Heat Release", "for " + ID]);

%% Calculate Heat of combustion and Temperature at exhaust - LHV way
% [T_exh, Q_combustion_LHV, m_combusted] = Calc_Q_LHV(C_p, mfr_fuel, mfr_air, RPM, W, LHV, T_int);

%% Calculate Heat of combustion and Temperature at exhaust - aROHR way
%[Q_combustion_aROHR] = Calc_Q_aROHR(p_filt, V_all, Ca, ValveEvents, gamma);
% disp(['Q combustion aROHR way: ', num2str(Q_combustion_aROHR), ' J']);

%% Thermal efficiency of the engine
% efficiency_LHV = (W / Q_combustion_LHV) *100; % efficiency for each cycle
% disp(['Calculated average thermal efficiency(LHV): ', num2str(efficiency_LHV), ' %']);
% 
% efficiency_LHV = (W / Q_combustion_aROHR) *100; % efficiency for each cycle
% disp(['Calculated average thermal efficiency(aROHR): ', num2str(efficiency_LHV), ' %']);

%% Calculate thermodynamic properties for each cycle - THIS STILL NEEDS TO BE IMPLEMENTED PROPERLY - NEED TO CALCULAT T_EXHAUST FOR IT SOMEHOW 
intake_species = [2, 3];           % Example species (O2 and N2)
exhaust_species = [4, 5, 3];       % Example species (CO2, H2O, and N2)
Y_int = [0.21, 0.79];              % Mole fractions for intake
Y_exh = [0.12, 0.18, 0.70];        % Mole fractions for exhaust

% % Call the function
%[Delta_H_all, Delta_U_all, Delta_S_all] = ThermoProperties(T_int, T_exh, SpS, Ncycles, Ca, intake_species, exhaust_species, Y_int, Y_exh);
 
% % Calculate averages
% Delta_H_avg = mean(Delta_H_all, 2);
% Delta_U_avg = mean(Delta_U_all, 2);
% Delta_S_avg = mean(Delta_S_all, 2);


%% Key performance indicators
% KPI data
% Format: data file, fuel, crank angle

% The selected fuel
MW_fuel = M_HVO;
% KPIdataFiles = HVO60_raw_dataFiles;
% Generate KPI table
KPITable = GenerateKPITable(IDsforKPI, true_mfr_fuel, T, LHV, RPM, AFR_stoich, x, MW_fuel,Cyl);
disp(KPITable)


%% Compare Raw and Filtered Pressure Data for a Single Cycle
figure;
hold on;
plot(Ca(:, 1), p_filt(:, 1) / bara, 'DisplayName', 'Raw Data');
plot(Ca(:, 1), p_avg(:, 1) / bara, 'DisplayName', 'Filtered Data');
xlabel('Crank Angle (°)');
ylabel('Pressure [bar]');
title(['Comparison of Raw and Filtered Pressure Data (Cycle ', num2str(iselect), ')']);
legend('show');
grid on;
hold off;

 


%% PV plot
[P_cycle, V_cycle] = IdealDieselCycle(Cyl) %Call ideal cycle function

figure; %Plot Ideal cycle
hold on;
plot(V_cycle, P_cycle / bara, 'DisplayName', 'ideal');
xlabel('Volume');
ylabel('Pressure [bar]');
title(['PV plot (Cycle ', num2str(iselect), ')']);
legend('show');
grid on;
hold off;


figure; %Plot experimental cycle
hold on;
plot(volume(:, 1), p_filt(:, 1) / bara, 'DisplayName', 'exp filtered data');
xlabel('Volume');
ylabel('Pressure [bar]');
title(['PV plot (Cycle ', num2str(iselect), ')']);
legend('show');
grid on;
hold off;
