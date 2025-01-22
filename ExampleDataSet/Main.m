% ===========================================
% MATLAB Script: Engine Data Analysis and Processing
% ===========================================
% This script loads experimental engine data and performs a series of analyses,
% including stoichiometric calculations, injection timing detection, heat release 
% computation, and plotting various performance graphs.
%
% Warning: All warnings are turned off in this script to avoid clutter.
% ===========================================

warning('off');  % Suppress warnings for cleaner output

%% Initialization
if ~exist("T","var")    % Ensure table T is defined
    load('T.mat','T')
end
clearvars -except T; close all; clc  % Clean workspace, close figures, clear command window
addpath("Functions", "Nasa"); savepath  % Add necessary paths for functions and NASA data

%% Units and Constants
% Conversion factors
mm = 1e-3;            % Millimeter to meter conversion
dm = 0.1;             % Decimeter to meter conversion
bara = 1e5;           % Bar to Pascals conversion
MJ = 1e6;             % MegaJoule to Joules conversion
kWhr = 1000 * 3600;   % Kilowatt-hour to Joules conversion
volperc = 0.01;       % Volume percentage conversion
ppm = 1e-6;           % Parts per million conversion

% Basic physical constants
g = 1e-3;             % Gram to kilogram conversion
s = 1;                % Time in seconds

% Experimental conditions
RPM = 1500;           % Constant RPM for experiments
T_int = 295.15;       % Assumed intake temperature [K]

% Molecular properties of fuels
M_diesel = 167;       % Molar mass of diesel [g/mol]
M_HVO = 226;          % Molar mass of HVO [g/mol]
M_GTL = 170;          % Molar mass of GTL [g/mol]

% Fuel densities [g/cm³]
Density_HVO = 0.78;
Density_Diesel = 0.85;
Density_GTL = 0.79;

% Lower heating values (LHV) [J/g]
LHV_diesel = 43e3;
LHV_HVO = 44.29e3;
LHV_GTL = 44e3;

% Carbon atom count in fuel molecules
x_HVO = 16;
x_diesel = 12;
x_GTL = 14;

%% Load and Reshape Data (also exhaust data)
ID = 'L50I17C50FuelHVO';  % Define experiment ID to load
% Load relevant data using custom function
[dataIn, ExhaustData, Ca, p_avg, S_current, mfr_fuel, CO_percent_load, HC_ppm_load, ...
 NOx_ppm_load, CO2_percent_load, O2_percent_load, lambda_load, x_blend, LHV_blend, ...
 perc_blend, fuel_used, M_blend] = loadingfromT(T, ID, bara, x_GTL, x_diesel, x_HVO, ...
                                               LHV_GTL, LHV_HVO, LHV_diesel, M_GTL, M_HVO, M_diesel);
IDsforKPI = {'L50I14C100FuelGTL', 'L50I16C100FuelGTL', 'L50I18C100FuelGTL', 'L50I20C100FuelGTL'};

%% Calculate LHV and Fuel Composition
perc_diesel = 1 - perc_blend;  
LHV = LHV_blend * perc_blend + LHV_diesel * perc_diesel;   % Fuel's lower heating value
x = x_blend * perc_blend + x_diesel * perc_diesel;       % Number of carbon atoms in fuel
MW_fuel = M_blend * perc_blend + perc_diesel * M_diesel; % Molar mass of fuel

%% Load NASA Data (if needed)
global Runiv
Runiv = 8.314;
[SpS, El] = myload('Nasa/NasaThermalDatabase.mat', {'N2', 'O2', 'CO2', 'H2O', 'Diesel'});

%% Compute Mass Flow
sps_fuel_name = 'Diesel';  % Fuel name for calculations
if mean(mfr_fuel) > 0.1 || mean(mfr_fuel) < 0.3
    mfr_fuel = mean(mfr_fuel);  % Use mean mass flow if within realistic bounds
    [~,~,AFR_stoich] = StoichiometricCombustion(SpS, El);
    mfr_air = CalculateMassFlowAir(O2_percent_load, mfr_fuel, AFR_stoich);
    AFR = mfr_air / mfr_fuel;
else
    % Use Bart's method if measured mfr_fuel is unrealistic
    M_C = 12; M_H = 1; M_CO2 = 44;
    mass_CH_ratio = 2;  
    CO2_mass_flow_rate = 0.5;
    molar_CH_ratio = mass_CH_ratio * (M_H / M_C);
    y = molar_CH_ratio * x;
    fuel_molar_mass = x * M_C + y * M_H;
    CO2_molar_flow_rate = CO2_mass_flow_rate / M_CO2;
    fuel_molar_flow_rate = CO2_molar_flow_rate / x;
    mfr_fuel = mean(fuel_molar_flow_rate * fuel_molar_mass);
end

%% Calculate Cylinder Volume
% Engine geometry parameters
Cyl.Bore = 104 * mm;
Cyl.Stroke = 85 * mm;
Cyl.CompressionRatio = 21.5;
Cyl.ConRod = 136.5 * mm;
Cyl.TDCangle = 180;
volume = CylinderVolume(Ca, Cyl);
disp('Cylinder volume calculated per cycle.');

%% Detect Injection Timing
% Compute average sensor current and set threshold
S_current_avg = mean(S_current, 2);
baseline_current = min(S_current_avg);
threshold = baseline_current + 0.2;

% Find injection start and end indices based on threshold
injection_indices = find(S_current_avg > threshold);
injection_start_idx = injection_indices(1);
injection_end_idx = injection_indices(end);

% Determine corresponding crank angles
injection_start_ca = Ca(injection_start_idx, 1);
injection_end_ca = Ca(injection_end_idx, 1);

% Display injection timing results
fprintf('Injection starts at %.2f° CA\n', injection_start_ca);
fprintf('Injection ends at %.2f° CA\n', injection_end_ca);

%% Plot Sensor Current with Injection Markers
figure;
plot(Ca(:, 1), S_current_avg, 'LineWidth', 1.5);
xlabel('Crank Angle (°)');
ylabel('Sensor Current (A)');
title('Average Sensor Current vs. Crank Angle');
grid on; xlim([-25, 25]);
hold on;
plot(injection_start_ca, S_current_avg(injection_start_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
plot(injection_end_ca, S_current_avg(injection_end_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
text(injection_start_ca, S_current_avg(injection_start_idx), sprintf('Start: %.2f°', injection_start_ca), ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
text(injection_end_ca, S_current_avg(injection_end_idx), sprintf('End: %.2f°', injection_end_ca), ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
hold off;

%% Plot Pressure vs. Crank Angle for All Cycles
figure;
plot(Ca, p_avg / bara, 'LineWidth', 1);
xlabel('Crank Angle (°)');
ylabel('Pressure (bar)');
xlim([-360, 360]); ylim([0, 50]);
title('Pressure vs. Crank Angle for All Cycles');
grid on;
iselect = 10;  % Specific cycle index to highlight
hold on;
plot(Ca(iselect), p_avg(iselect) / bara, 'r', 'LineWidth', 2);
set(gca, 'XTick', -360:60:360);
grid on;

%% Calculate Work
W = trapz(volume, p_avg * 1e5);  % Compute work from pressure-volume curve
disp(['Calculated work: ', num2str(W), ' J']);

%% Stoichiometric Calculations for Diesel and HVO
% Stoichiometric combustion for diesel
[stoich_coeffs, reaction_eq, AFR_stoich] = StoichiometricCombustion(SpS, El);
stoich_coeffs.fuel

% Stoichiometric calculations for HVO
m_fuel_per100cycles = 0.128;  % Fuel mass per 100 cycles [g]
moles_HVO = m_fuel_per100cycles / M_HVO;
moles_O2 = moles_HVO * (3*x_HVO*0.5) + 1/2;
moles_N2 = moles_O2 * 79/21;
moles_CO2 = moles_HVO * 16;
moles_H2O = moles_HVO * 17;
Stoich_HVO100 = [moles_HVO, moles_O2, moles_N2, moles_CO2, moles_H2O, moles_N2];

% Calculations for HVO50 blend
a = m_fuel_per100cycles / (((Density_HVO / Density_Diesel) + 1) * M_diesel);
b = m_fuel_per100cycles / (((Density_Diesel / Density_HVO) + 1) * M_HVO);
e = a * 12 + b * 16;
g = a * 23 + b * 34;
q = 9.25 * a + 24.5 * b;
k = q * (79 / 21);
Stoich_HVO50 = [a, b, q, k, e, g/2, k];

%% Compute Apparent Rate of Heat Release (aROHR) and Heat Release (aHR)
if exist('O2_percent_load','var')
    gamma = CalculateGamma(SpS, volume, p_avg, O2_percent_load, CO2_percent_load, mfr_fuel, AFR_stoich, RPM);
else
    gamma = 1.32;  % Default value if exhaust data is unavailable
end
aROHR = get_aROHR(p_avg, volume, gamma);
% Isolate peak region for analysis
idxStart = 355 / 0.2; 
idxEnd = (50 + 360) / 0.2;
aROHR(1:idxStart) = 0; 
aROHR(idxEnd:end) = 0;
aHR = get_aHR(aROHR);
disp("aROHR and aHR computed successfully!");

%% Plot Apparent Heat Release
figure;
subplot(1,2,1)
plot(Ca, aROHR);
xlabel("Crank Angle [deg]");
ylabel("Apparent Rate of Heat Release [J/deg]");
xlim([-10,30]);
legend("aROHR","Location","southeast");
title(["Apparent Rate of Heat", "Release for " + fuel_used]);

subplot(1,2,2)
plot(Ca, aHR);
xlabel("Crank Angle [deg]");
ylabel("Apparent Heat Release [J]");
xlim([-10,30]);
legend("aHR","Location","southeast");
title(["Apparent Heat Release", "for " + fuel_used]);

figure;
subplot(1,2,1)
plot(Ca, aROHR);
xlabel("Crank Angle [deg]");
ylabel("Apparent Rate of Heat Release [J/deg]");
xlim([-5,50]);
legend("aROHR","Location","southeast");
title(["Apparent Rate of Heat", "Release for " + fuel_used]);

subplot(1,2,2)
plot(Ca, aHR);
xlabel("Crank Angle [deg]");
ylabel("Apparent Heat Release [J]");
xlim([-5,50]);
legend("aHR","Location","southeast");
title(["Apparent Heat Release", "for " + fuel_used]);

%% Key Performance Indicators
% Generate KPI table using custom function
KPITable = GenerateKPITable(IDsforKPI, mfr_fuel, T, LHV, RPM, AFR_stoich, x, MW_fuel, Cyl, fuel_used);
disp(KPITable);

%% Compare Raw and Filtered Pressure Data for a Single Cycle
figure;
hold on;
plot(Ca(:, 1), p_avg(:, 1) / bara, 'DisplayName', 'Raw Data');
plot(Ca(:, 1), p_avg(:, 1) / bara, 'DisplayName', 'Filtered Data');
xlabel('Crank Angle (°)');
ylabel('Pressure [bar]');
title(['Comparison of Raw and Filtered Pressure Data (Cycle ', num2str(iselect), ')']);
legend('show');
grid on;
hold off;

%% PV Plot: Experimental vs. Ideal Cycle
p_avg = (p_avg * 10e4);  % Convert pressure units
[P_cycle, V_cycle] = IdealDieselCycle(Cyl);  % Get ideal cycle data

figure;
hold on;
plot(V_cycle, P_cycle / bara, 'DisplayName', 'Ideal Cycle');
plot(volume(:, 1), p_avg(:, 1) / bara, 'DisplayName', 'Experimental Data');
xlabel('Volume');
ylabel('Pressure [bar]');
title(['PV Plot (Cycle ', num2str(iselect), ')']);
legend('show');
grid on;
hold off;

figure;
hold on;
plot(V_cycle, P_cycle / bara, 'DisplayName', 'Ideal Cycle');
plot(volume(:, 1), p_avg(:, 1) / bara, 'DisplayName', 'Experimental Data');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Volume');
ylabel('Pressure [bar]');
title(['PV Plot (Cycle ', num2str(iselect), ') - Log Scale']);
legend('show');
grid on;
hold off;

% ===========================================
% End of Script
% ===========================================
