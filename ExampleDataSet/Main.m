warning off
%% Initialization
clear; clc; close all;
%addpath("Functions", "Nasa");  % Add necessary paths
figure('Visible', 'on');
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
RPM = 1500;     % constant RPM of experiments
M_diesel = 167; % Molar mass (g/mol)
M_HVO = 226; % Molar mass (g/mol)
Density_HVO = 0.78;
Density_Diesel = 0.85;
%% Making Path for files of different blends

% Format: data file, fuel, crank angle
HVO50_raw_dataFiles = {
        fullfile('Data', 'session2_Raw', '20241212_0000002_HVO50_secondexperiment_CA14.txt'), 'HVO50', 14;
        fullfile('Data', 'session2_Raw', '20241212_0000008_HVO50_secondexperiment_CA15.txt'), 'HVO50', 15;
        fullfile('Data', 'session2_Raw', '20241212_0000003_HVO50_secondexperiment_CA16.txt'), 'HVO50', 16;
        fullfile('Data', 'session2_Raw', '20241212_0000004_HVO50_secondexperiment_CA17.txt'), 'HVO50', 17;
        fullfile('Data', 'session2_Raw', '20241212_0000005_HVO50_secondexperiment_CA18.txt'), 'HVO50', 18;
        fullfile('Data', 'session2_Raw', '20241212_0000009_HVO50_secondexperiment_CA19.txt'), 'HVO50', 19;
        fullfile('Data', 'session2_Raw', '20241212_0000006_HVO50_secondexperiment_CA20.txt'), 'HVO50', 20;

    };

HVO60_raw_dataFiles = {
        fullfile('Data', 'Load50-MultipleSOI', 'SOI14.txt'), 'HVO60', 14;
        fullfile('Data', 'Load50-MultipleSOI', 'SOI15.txt'), 'HVO60', 15;
        fullfile('Data', 'Load50-MultipleSOI', 'SOI16.txt'), 'HVO60', 16;
        fullfile('Data', 'Load50-MultipleSOI', 'SOI17.txt'), 'HVO60', 17;
        fullfile('Data', 'Load50-MultipleSOI', 'SOI18.txt'), 'HVO60', 18;
        fullfile('Data', 'Load50-MultipleSOI', 'SOI19.txt'), 'HVO60', 19;
        fullfile('Data', 'Load50-MultipleSOI', 'SOI20.txt'), 'HVO60', 20;

    };

Processed_session1_files = {
    fullfile('Data', 'processed_session1', '15CA_filtered_averaged.txt'), 'Diesel', 15;
    fullfile('Data', 'processed_session1', '16CA_filtered_averaged.txt'), 'Diesel', 16;
    fullfile('Data', 'processed_session1', '17CA_filtered_averaged.txt'), 'Diesel', 17;
    fullfile('Data', 'processed_session1', '18CA_filtered_averaged.txt'), 'Diesel', 18;
    fullfile('Data', 'processed_session1', '19CA_filtered_averaged.txt'), 'Diesel', 19;
    fullfile('Data', 'processed_session1', '20CA_filtered_averaged.txt'), 'Diesel', 20;
    fullfile('Data', 'processed_session1', '21CA_filtered_averaged.txt'), 'Diesel', 21;
};

session1_Raw_files = {
    fullfile('Data', 'session1_Raw', '20241125_0000010_15CA.txt'), 'Diesel', 15;
    fullfile('Data', 'session1_Raw', '20241125_0000014_16CA.txt'), 'Diesel', 16;
    fullfile('Data', 'session1_Raw', '20241125_0000016_17CA.txt'), 'Diesel', 17;
    fullfile('Data', 'session1_Raw', '20241125_0000013_18CA.txt'), 'Diesel', 18;
    fullfile('Data', 'session1_Raw', '20241125_0000011_19CA.txt'), 'Diesel', 19;
    fullfile('Data', 'session1_Raw', '20241125_0000012_20CA.txt'), 'Diesel', 20;
    fullfile('Data', 'session1_Raw', '20241125_0000018_21CA.txt'), 'Diesel', 21;
};

%% Define Fuel used
fuel_name = 'Diesel';
LHV = 43e3; %Lower heating value given in the project guide for Diesel B7 J/g
O2_perc = 14.42; % O2 percentage at exhaust (hardcoded)
x = 12;


%% Load NASA Data (if needed)
global Runiv
Runiv = 8.314;

[SpS, El] = myload('Nasa/NasaThermalDatabase.mat', {'Diesel', 'O2', 'N2', 'CO2', 'H2O'});

%% Engine Geometry Data
Cyl.Bore = 104 * mm;
Cyl.Stroke = 85 * mm;
Cyl.CompressionRatio = 21.5;
Cyl.ConRod = 136.5 * mm;
Cyl.TDCangle = 180;

%% Valve Events (Crank Angles)
ValveEvents.CaIVO = -355;
ValveEvents.CaIVC = -135;
ValveEvents.CaEVO = 149;
ValveEvents.CaEVC = -344;
ValveEvents.CaSOI = -3.2;  % Start of Injection

%% Load and Reshape Data
dataFileName = fullfile('Data' , 'processed_Data_experiment1_load3.5.txt');
dataIn = table2array(readtable(dataFileName));

%Error handling:
[Nrows, Ncols] = size(dataIn);
if Ncols  ~= 4
    warning('Data loaded does not have the expected amount of columns (4)');
end

resolution = 0.2;  % Degrees crank angle resolution
NdatapointsPerCycle = 720 / resolution;
Ncycles = Nrows / NdatapointsPerCycle;

if mod(Nrows, NdatapointsPerCycle) ~= 0
    error('Number of data points is not an integer multiple of data points per cycle.');
end

% Reshape data into cycles
Ca = reshape(dataIn(:, 1), [], Ncycles);      % Crank angle in degrees
p = reshape(dataIn(:, 2), [], Ncycles) * bara;  % Pressure in Pa
S_current = reshape(dataIn(:, 3), [], Ncycles);  % Sensor current 
mfr_fuel = reshape(dataIn(:, 4), [], Ncycles);  % Fuel mass flow


%% Filter Pressure Data
polynomialOrder = 3;
frameLength = 21;  % Must be odd

% Initialize the filtered pressure matrix
p_filtered = zeros(size(p));

% Apply the filter to each column
for i = 1:Ncycles
    p_filtered(:, i) = SGFilter(p(:, i), polynomialOrder, frameLength, 0);
end

disp('Data filtered and reshaped into cycles');

%% load the excelfile
fileName = fullfile('Data', 'Session2.xlsx');
sheetName = 'Sheet2';
range1 = 'A8:G16'; % Table 1 range
range2 = 'A22:G28'; % Table 2 range

% Read Table 1
table1experiment1 = readtable(fileName, 'Sheet', sheetName, 'Range', range1);

load_data = table1experiment1{:, 1};      % Load
CO_percent_load = table1experiment1{:, 2}; % CO%
HC_ppm_load = table1experiment1{:, 3};    % HC ppm
NOx_ppm_load = table1experiment1{:, 4};   % NOx ppm
CO2_percent_load = table1experiment1{:, 5}; % CO2%
O2_percent_load = table1experiment1{:, 6}; % O2%
lambda_load = table1experiment1{:, 7};    % Lambda

% Read Table 2
table2experiment1 = readtable(fileName, 'Sheet', sheetName, 'Range', range2);
table2experiment1.Properties.VariableNames = {'CrankAngle', 'CO_percent', 'HC_ppm', 'NOx_ppm', 'CO2_percent', 'O2_percent', 'Lambda'};
disp('Table 2:');
disp(table2experiment1);

CA = table2experiment1{:, 1};            % Crank Angle (CA)
CO_percent_CA = table2experiment1{:, 2}; % CO%
HC_ppm_CA = table2experiment1{:, 3};    % HC ppm
NOx_ppm_CA = table2experiment1{:, 4};   % NOx ppm
CO2_percent_CA = table2experiment1{:, 5}; % CO2%
O2_percent_CA = table2experiment1{:, 6}; % O2%
lambda_CA = table2experiment1{:, 7};    % Lambda


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
set(gcf, 'Position', [200, 800, 1200, 400]);

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
set(gcf, 'Position', [200, 800, 1200, 400]);

plot(Ca, p_filtered / bara, 'LineWidth', 1);
xlabel('Crank Angle (°)');
ylabel('Pressure (bar)');
xlim([-360, 360]);
ylim([0, 50]);
title('Pressure vs. Crank Angle for All Cycles');
grid on;

% Highlight a specific cycle
iselect = 10;
hold on;
plot(Ca(:, iselect), p(:, iselect) / bara, 'r', 'LineWidth', 2);

% Plot valve events
YLIM = ylim;
line([ValveEvents.CaIVC, ValveEvents.CaIVC], YLIM, 'Color', 'b', 'LineWidth', 1);
line([ValveEvents.CaEVO, ValveEvents.CaEVO], YLIM, 'Color', 'r', 'LineWidth', 1);

set(gca, 'XTick', -360:60:360);
grid on;

%% Calculate Cylinder Volume for All Cycles
V_all = zeros(size(Ca));  % Initialize volume matrix
for i = 1:Ncycles
    V_all(:, i) = CylinderVolume(Ca(:, i), Cyl);
end

disp('Cylider volume calculated / cycle');

%% Calculate Average Volume and Pressure
V_avg = mean(V_all, 2);         % Average volume across all cycles for every CA
p_avg = mean(p, 2);             % Average pressure across all cycles for every CA
p_filtered_avg = mean(p_filtered, 2);

%% Calculate Work
W = trapz(V_avg, p_avg); % Calculate the area under the averaged p-V curve
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

%% Stoichiometric calculations
[stoich_coeffs, reaction_eq, AFR_stoich] = StoichiometricCombustion(fuel_name, SpS, El);
stoich_coeffs.fuel
%% Calculate mass flow of air:
O2_percent = mean(O2_percent_load(4:6));  % Load relevant exhaust data from processed excel
mfr_air = CalculateMassFlowAir(O2_percent, true_mfr_fuel, AFR_stoich);
fprintf('Mass flow rate for air: %.6f g/s\n', mfr_air)
%% Calculate Cp and gamma
% [cp, gamma] = calc_cp_gamma(LHV, mfr_fuel, mfr_air);

%% Cp and gamma
C_p =  1101.6; %Cp manually plugged in from the results of cp and gamma calculations [J/g*K]
gamma = 1.312562;  % gamma manually plugged in from the results of cp and gamma calculations

%% Calculate Temperature at exhaust
T_int = 295.15 * ones(1, 100); %assume ambient intake temperature (22C) [K]
[rowsmfr_fuel, colsmfr_fuel] = size(mfr_fuel);

[T_exh, Q_combustion_percycle, avg_m_fuelpercycle] = Texhaust(CA, C_p, mfr_fuel, mfr_air, RPM, W, LHV, T_int, Ncycles);

%% Run the energy of combustion calculation - aROHR way
[Q_combustion] = MASSHeatOfCombustion(p_filtered, V_all, Ca, ValveEvents, gamma);
disp(['Q combustion aROHR way: ', num2str(Q_combustion), ' J']);

%% Thermal efficiency of the engine
efficiency_percycle = W ./ Q_combustion_percycle; % efficiency for each cycle
efficiency_avg = mean(efficiency_percycle(2:end)) *100; % average efficieny of cycles (the first cycle is excluded as it varies from the rest - due to starting up)

disp(['Calculated average thermal efficiency: ', num2str(efficiency_avg), ' %']);

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
plot(Ca(:, iselect), p_filtered(:, iselect) / bara, 'DisplayName', 'Filtered Data');
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
 