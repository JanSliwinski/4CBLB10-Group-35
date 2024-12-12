warning off
%% Initialization
clear; clc; close all;
%addpath("Functions", "Nasa");  % Add necessary paths

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

%% Define Fuel used
fuel_name = 'Diesel';
LHV = 43e3; %Lower heating value given in the project guide for Diesel B7 J/g
x = 12;
MW_Fuel = 200; %Molar weight of fuel


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

%% Process experimental data
folderpath = './ExampleDataSet/Data/session1_Raw/load3.5'; % path to raw data
outputfilePath = './ExampleDataSet/Data/processed_Data_experiment1_load3.5.txt'; % path to output file
% averagedata = AverageExperimentData(folderpath, outputfilePath); % run function averaging relevant data

%% Load and Reshape Data
dataFileName = fullfile('Data' ,'processed_Data_experiment1_load3.5.txt');
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
% mfr_fuel = reshape(dataIn(:, 4), [], Ncycles);  % Fuel mass flow
mfr_fuel = 0.16;    % assumed constant value for mass flow of fuel [g/s]

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
fileName = fullfile('Data', 'Session1.xlsx');
sheetName = 'Sheet1';
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
disp('Table 2:');
disp(table2experiment1);

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

%% Calculate mass flow of air:
O2_percent = mean(O2_percent_load(4:6));  % Load relevant exhaust data from processed excel
mfr_air = CalculateMassFlowAir(O2_percent, mfr_fuel, AFR_stoich);

%% Calculate Cp and gamma
% [cp, gamma] = calc_cp_gamma(LHV, mfr_fuel, mfr_air);

%% Cp and gamma
%C_p = 1005; % Cp of air
C_p =  1101.6; %Cp manually plugged in from the results of cp and gamma calculations [J/g*K]
gamma = 1.312562;  % gamma manually plugged in from the results of cp and gamma calculations

%% Calculate Heat of combustion and Temperature at exhaust - LHV way
[T_exh, Q_combustion_LHV, m_combusted] = Calc_Q_LHV(C_p, mfr_fuel, mfr_air, RPM, W, LHV, T_int);

%% Calculate Heat of combustion and Temperature at exhaust - aROHR way
[Q_combustion_aROHR] = Calc_Q_aROHR(p_filtered, V_all, Ca, ValveEvents, gamma);
disp(['Q combustion aROHR way: ', num2str(Q_combustion_aROHR), ' J']);

%% Thermal efficiency of the engine
efficiency_LHV = (W / Q_combustion_LHV) *100; % efficiency for each cycle
disp(['Calculated average thermal efficiency(LHV): ', num2str(efficiency_LHV), ' %']);
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


KPIdataFiles = {
        fullfile('Data', 'session1_Raw', '20241125_0000010_15CA.txt'), 'Diesel', 15;
        fullfile('Data', 'session1_Raw', '20241125_0000014_16CA.txt'), 'Diesel', 16;
        fullfile('Data', 'session1_Raw', '20241125_0000016_17CA.txt'), 'Diesel', 17;
        fullfile('Data', 'session1_Raw', '20241125_0000013_18CA.txt'), 'Diesel', 18;
        fullfile('Data', 'session1_Raw', '20241125_0000011_19CA.txt'), 'Diesel', 19;
        fullfile('Data', 'session1_Raw', '20241125_0000012_20CA.txt'), 'Diesel', 20;
        fullfile('Data', 'session1_Raw', '20241125_0000017_21CA.txt'), 'Diesel', 21;
    };

 %Generate KPI table
KPITable = GenerateKPITable(KPIdataFiles, table2experiment1, LHV, avg_m_fuelpercycle, RPM, AFR_stoich, x, MW_Fuel,Cyl);
disp(KPITable)

%% Rate of changes, Pressure and Volume
% Crank angle change per data point
dCA = resolution;

% Pressure change per data point
dp = diff(p_filtered_avg);
% Pressure change per Crank angle
dp_dCA = dp/dCA;

% Volume change per data point
dV = diff(V_avg);
% Volume change per Crank angle
dV_dCA = dV/dCA;

%% Calculate aROHR
 
aROHR_avg = aROHR(p_filtered_avg, V_avg, resolution, gamma, dp_dCA, dV_dCA);

% Plot the apparent Rate of Heat Release
figure;
plot(Ca(:, 1), aROHR_avg, 'LineWidth', 1.5);
xlabel('Crank Angle (°)');
ylabel('Apparent Rate of Heat Release [J/°]');
title('Apparent Rate of Heat Release (Average)');
xlim([-35,135])
grid on;

%% Calculate Apparent Heat Release
aHR_avg = aHR(aROHR_avg, resolution);  % Assuming aHR function is already defined

% Plot the Apparent Heat Release
figure;
plot(Ca(:, 1), aHR_avg, 'LineWidth', 1.5);  % Ca(:,1) is the crank angle array
xlabel('Crank Angle (°)');
ylabel('Apparent Heat Release [J]');
title('Apparent Heat Release (Average)');
xlim([-45,135]);
grid on;
hold on;

% Determine the Maximum Value and Its Index
[max_aHR, max_idx] = max(aHR_avg);  % Find the max value and its index

% Calculate 10%, 50%, and 90% of Maximum Value
val_10 = 0.1 * max_aHR;  % 10% of max
val_50 = 0.5 * max_aHR;  % 50% of max
val_90 = 0.9 * max_aHR;  % 90% of max

% Find Crank Angles Before the Peak
% Use only the range before and including the peak for interpolation
crank_angle_pre_peak = Ca(1:max_idx, 1);  % Crank angles up to the peak
aHR_pre_peak = aHR_avg(1:max_idx);       % aHR values up to the peak

% Interpolate for the 10%, 50%, and 90% values
crank_angle_10 = interp1(aHR_pre_peak, crank_angle_pre_peak, val_10);
crank_angle_50 = interp1(aHR_pre_peak, crank_angle_pre_peak, val_50);
crank_angle_90 = interp1(aHR_pre_peak, crank_angle_pre_peak, val_90);

% Plot the Updated Results
% Highlight points with scatter
scatter([crank_angle_10, crank_angle_50, crank_angle_90], ...
        [val_10, val_50, val_90], 'r', 'filled');

% Annotate the points
text(crank_angle_10, val_10, sprintf('10%% (%.2f°)', crank_angle_10), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
text(crank_angle_50, val_50, sprintf('50%% (%.2f°)', crank_angle_50), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
text(crank_angle_90, val_90, sprintf('90%% (%.2f°)', crank_angle_90), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);


%% Plot pV Diagrams
figure;
tl = tiledlayout(2, 2);

% Subplot 1: Linear pV Diagram (Average)
nexttile;
plot(V_avg / dm^3, p_avg / bara);
xlabel('Volume [dm³]');
ylabel('Pressure [bar]');
title('Average pV Diagram (Linear Scale)');
grid on;

% Subplot 2: Logarithmic pV Diagram (Average)
nexttile;
loglog(V_avg / dm^3, p_avg / bara);
xlabel('Volume [dm³]');
ylabel('Pressure [bar]');
title('Average pV Diagram (Logarithmic Scale)');
grid on;

% Subplot 3: Single Cycle pV Diagram
nexttile;
plot(V_all(:, iselect) / dm^3, p(:, iselect) / bara);
xlabel('Volume [dm³]');
ylabel('Pressure [bar]');
title(['pV Diagram (Cycle ', num2str(iselect), ')']);
grid on;

% Subplot 4: Filtered Average pV Diagram
nexttile;
loglog(V_avg / dm^3, p_filtered_avg / bara);
xlabel('Volume [dm³]');
ylabel('Pressure [bar]');
title('Filtered Average pV Diagram (Logarithmic Scale)');
grid on;

sgtitle('pV Diagrams');

if any(p(:) < 0)
    warning('There are negative pressure values.');
else
    disp('All pressure values are non-negative.');
end

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

%% Ideal Diesel Cycle
% Extract initial conditions from actual data
% Find index corresponding to Intake Valve Closure (IVC)
% IVC_angle = ValveEvents.CaIVC;  % Crank angle for IVC
% P1 = p_filtered_avg(idx_IVC);  % Pressure at IVC
% V1 = V_avg(idx_IVC);           % Volume at IVC
% Assume inlet temperature or estimate based on conditions
% <<<<<<< HEAD:ExampleDataSet/Simple.m
% T1 = 295 ;  % K assumed to be ambient
% =======
% T1 = 298.15;  % K (Atmospheric Conditions)




