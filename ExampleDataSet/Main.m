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
% LHV_HVO = ?

%% Load and Reshape data (also exhaust data)
ID = 'L50I17C50'; %DEFINE ID OF THE EXPERIMENT DATA YOU WANT TO LOAD IN!
%run fucntion to load in all relevant data

[dataIn, ExhaustData, Ca, p_filt, S_current, mfr_fuel, CO_percent_load, HC_ppm_load, NOx_ppm_load, CO2_percent_load, O2_percent_load, lambda_load] = loadingfromT(T, ID, bara);


%% Define Fuel used
fuel_name = 'Diesel';
fprintf('Fuel used:', fuel_name);
LHV = LHV_diesel;
x = 12;     %Whats this???? 

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

%% Calculate Work
W = trapz(volume, p_filt*1e5); % Calculate the area under the averaged p-V curve
disp(['Calculated work: ', num2str(W), ' J']);

% Constants
M_C = 12; % Molar mass of Carbon (g/mol)
M_H = 1; % Molar mass of Hydrogen (g/mol)
M_CO2 = 44; % Molar mass of CO2 (g/mol)
mass_CH_ratio = 5.49; % Typical C/H mass ratio for diesel

% Inputs: Known CO2 mass flow rate
CO2_mass_flow_rate = 0.5; % IDK what this value is but Barrt Sommers says we should have it

% Convert mass CH ratio to molar CH ratio
molar_CH_ratio = mass_CH_ratio * (M_H / M_C);

% Determine x and y for diesel
x = 12; % Carbon atoms in standard diesel
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

%% Calculate mass flow of air:

% mfr_air = CalculateMassFlowAir(O2_percent, mfr_fuel, AFR_stoich);
% mfr_air = CalculateMassFlowAir(O2_percent, true_mfr_fuel, AFR_stoich);
% fprintf('Mass flow rate for air: %.6f g/s\n', mfr_air)


%% aROHR
p_filtRough = sgolayfilt(p_filt,2,101);
if exist('O2_percent_load','var')
    gamma = CalculateGamma(SpS,volume,p_filtRough,O2_percent_load,CO2_percent_load,true_mfr_fuel,AFR_stoich,RPM);
else
    gamma = 1.32; % constant if no exhuast data exist
end
aROHR = get_aROHR(p_filtRough,volume,gamma);
% Isolate peak
idxStart = 355 / 0.2; idxEnd = (50 + 360) / 0.2;
aROHR(1:idxStart) = 0; aROHR(idxEnd:end) = 0;
aHR = get_aHR(aROHR);
disp("aROHR and aHR computed successfully!")

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
KPIdataFiles = HVO60_raw_dataFiles;
% Generate KPI table
KPITable = GenerateKPITable(KPIdataFiles, true_mfr_fuel, table2experiment1, LHV, RPM, AFR_stoich, x, MW_fuel,Cyl);
disp(KPITable)


%% Compare Raw and Filtered Pressure Data for a Single Cycle
figure;
hold on;
plot(Ca(:, iselect), p(:, iselect) / bara, 'DisplayName', 'Raw Data');
plot(Ca(:, iselect), p_filt(:, iselect) / bara, 'DisplayName', 'Filtered Data');
xlabel('Crank Angle (°)');
ylabel('Pressure [bar]');
title(['Comparison of Raw and Filtered Pressure Data (Cycle ', num2str(iselect), ')']);
legend('show');
grid on;
hold off;
%% Stoichiometric Calculations for HVO
x_HVO = 16;
m_fuel_per100cycles = 0.128; %(g)
%HVO100
moles_HVO = m_fuel_per100cycles / M_HVO;

moles_O2 = moles_HVO * (3*x_HVO*0.5) + 1/2; % Moles of O2 required
moles_N2 = moles_O2 * 79/21;   % Moles of N2 in reactants


moles_CO2 = moles_HVO * 16; % Moles of CO2 produced
moles_H2O = moles_HVO * 17; % Moles of H2O produced
Stoich_HVO100 = [moles_HVO, moles_O2, moles_N2, moles_CO2, moles_H2O, moles_N2];
%HVO50
m_fuel_per100cycles = 0.128; %(g)

 a = m_fuel_per100cycles / (((Density_HVO / Density_Diesel) + 1) * M_diesel);
 b = m_fuel_per100cycles / (((Density_Diesel / Density_HVO) + 1) * M_HVO);
 e = a * 12 + b * 16;
 g = a * 23 + b * 34;
 q = 9.25 * a + 24.5 * b;
 k = q * (79 / 21);
 Stoich_HVO50 = [a, b, q, k, e, g/2, k]; %Moles of each Component per 100 cycles
 