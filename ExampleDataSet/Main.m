warning off
%% Initialization
clear; clc; close all;
addpath("Functions", "Nasa");  % Add necessary paths

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
%ProcessExperimentData('./session1_Raw/Group 35', 'processed_Data.txt');

%% Load and Reshape Data
dataFileName = fullfile('Data' ,'processed_Data_experiment1.txt');
dataIn = table2array(readtable(dataFileName));

[Nrows, ~] = size(dataIn);
resolution = 0.2;  % Degrees crank angle resolution
NdatapointsPerCycle = 720 / resolution;
Ncycles = Nrows / NdatapointsPerCycle;

if mod(Nrows, NdatapointsPerCycle) ~= 0
    error('Number of data points is not an integer multiple of data points per cycle.');
end

% Reshape data into cycles
Ca = reshape(dataIn(:, 1), [], Ncycles);      % Crank angle in degrees
p = reshape(dataIn(:, 2), [], Ncycles) * bara;  % Pressure in Pa
SC = reshape(dataIn(:, 3), [], Ncycles);  % Sensor current [mA] 
mfr_fuel = reshape(dataIn(:, 4), [], Ncycles);  % Fuel mass flow

%% load the excelfile
fileName = 'Session1.xlsx';
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

plot(Ca, p / bara, 'LineWidth', 1);
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

%% Filter Pressure Data
polynomialOrder = 3;
frameLength = 21;  % Must be odd

% Initialize the filtered pressure matrix
p_filtered = zeros(size(p));

% Apply the filter to each column
for i = 1:Ncycles
    p_filtered(:, i) = SGFilter(p(:, i), polynomialOrder, frameLength, 0);
end

%% Calculate Average Volume and Pressure
V_avg = mean(V_all, 2);         % Average volume across all cycles
p_avg = mean(p, 2);             % Average pressure across all cycles
p_filtered_avg = mean(p_filtered, 2);

%% Calculate Work
W = trapz(V_avg, p_avg);
disp(['Calculated work: ', num2str(W), ' J']);

%% Calculate thermodynamic properties for each cycle
% intake_species = [2, 3];           % Example species (O2 and N2)
% exhaust_species = [4, 5, 3];       % Example species (CO2, H2O, and N2)
% Y_int = [0.21, 0.79];              % Mole fractions for intake
% Y_exh = [0.12, 0.18, 0.70];        % Mole fractions for exhaust
% 
% % Call the function
% [Delta_H_all, Delta_U_all, Delta_S_all] = ThermoProperties(T_int, T_exh, SpS, Ncycles, Ca, intake_species, exhaust_species, Y_int, Y_exh);
% 
% % Calculate averages
% Delta_H_avg = mean(Delta_H_all, 2);
% Delta_U_avg = mean(Delta_U_all, 2);
% Delta_S_avg = mean(Delta_S_all, 2);

%% Key performance indicators
% % Power calculation
% P = W*(RPM/2*60);
% 
% % Calls the KPI function
% KPIs = CalculateKPIs(W, mfr_fuel, LHV, P, mfr_CO2, mfr_NOx);


%% Load Exhaust Gas Data

O2_percent_vector = repmat(O2_percent_load, 1, 100);

%% Calculate Air Mass Flow Rate
AFR_stoich = 14.5;  % Stoichiometric AFR for diesel
mfr_air = CalculateMassFlowAir(O2_percent_vector, mfr_fuel, AFR_stoich);

disp('Air mass flow rate calculated successfully.');



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



