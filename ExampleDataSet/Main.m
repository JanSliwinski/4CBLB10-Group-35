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
RPM = 1500;     %constant RPM of experiments

%% Define Fuel used
fuel_name = 'Diesel';
LHV = 43e3; %Lower heating value given in the project guide for Diesel B7 J/g
O2_perc = 14.42; % O2 percentage at exhaust (hardcoded)

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
averagedata = AverageExperimentData(folderpath, outputfilePath); % run function averaging relevant data

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
%<<<<<<< HEAD
disp('Table 2:');
disp(table2experiment1);
%>>>>>>> b2452164e72c33b9a41ecc6c9454d555500da4ba
%=======

CA = table2experiment1{:, 1};            % Crank Angle (CA)
CO_percent_CA = table2experiment1{:, 2}; % CO%
HC_ppm_CA = table2experiment1{:, 3};    % HC ppm
NOx_ppm_CA = table2experiment1{:, 4};   % NOx ppm
CO2_percent_CA = table2experiment1{:, 5}; % CO2%
O2_percent_CA = table2experiment1{:, 6}; % O2%
lambda_CA = table2experiment1{:, 7};    % Lambda
%>>>>>>> 51126a27191a1963030fa9c8e65cc87d3e1a906e

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

%% Stoichiometric calculations
[stoich_coeffs, reaction_eq, AFR_stoich] = StoichiometricCombustion(fuel_name, SpS, El);

%% Calculate mass flow of air:
O2_percent = O2_perc;
mfr_air = CalculateMassFlowAir(O2_percent, mfr_fuel, AFR_stoich);

%% Calculate Temperature at exhaust
T_int = 295.15 * ones(1, 100); %assume ambient intake temperature (22C) [K]
[rowsmfr_fuel, colsmfr_fuel] = size(mfr_fuel);

% Calculate time step per crank angle
DeltaCA = mean(diff(CA(:,1))); % Step size between crank angles (degrees)
Delta_t_perCA = (60 * DeltaCA) / (RPM * 360); % Time per crank angle step (s)
m_fuel_percycle = sum(mfr_fuel, 1) * Delta_t_perCA; % calculates the total mass of fuel in /cycle [g]

Q_combustion_percycle = m_fuel_percycle * LHV; % energy of combustion /cycle [J]
% But some part of this energy will be used as work
Q_warmingtheair = Q_combustion_percycle - W; % energy used for warming the exhaust gasses [J]

% The rest is used for warming the exhaust gasses
mfr_exh = mfr_fuel + mfr_air;    % Mass flow rate at the exhaust [g/s]
[rowsmfr_exh, colsmfr_exh] = size(mfr_exh);
m_exh_percycle = sum(mfr_exh, 1) * Delta_t_perCA; % the total mass at exhaust /cycle [g]
C_p = 1.005; %assume Cp of air [J/g*K]
delta_T_percycle = Q_warmingtheair ./ (C_p * m_exh_percycle); % change in temperature due to combustion/cycle    
T_exh = T_int + delta_T_percycle; % exhaust temperature /cycle


%checking whether data is realistic
avg_m_fuelpercycle = mean(m_fuel_percycle(2:end));
avg_m_exh_percycle = mean(m_exh_percycle(2:end));
avg_delt_T = mean(delta_T_percycle(2:end));
avg_T_exh = mean(T_exh(2:end));
avg_Q = mean(Q_combustion_percycle(2:end));
disp(['Average fuel mass per cycle: ', num2str(avg_m_fuelpercycle)]);
%disp(['Average exhaust mass per cycle: ', num2str(avg_m_exh_percycle)]);
disp(['Average energy of combustion: ', num2str(avg_Q)]);
disp(['Average delta T: ', num2str(avg_delt_T)]);
disp(['Average temperature at exhaust: ', num2str(avg_T_exh)]);

% Jings code for efficiency
efficiency = (W / (mean(mfr_fuel(:)) * LHV))*100;
disp(['Jing-s efficiency: ', num2str(efficiency), ' %']);

% Plot the change of T_exhaust over all cycles
figure;
cycles = 1:Ncycles;
plot(cycles, T_exh, 'LineWidth', 1);
xlabel('Cycles');
ylabel('Exhaust temperature (K)');
xlim([1, 100]);
title('Exhaust temperature over the 100 cycles run');
grid on;

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

%% Key performance indicators - this still needs to be implemented properly!!
%mfr_CO2 = 
%mfr_Nox = 

% % Power calculation
%<<<<<<< HEAD
P = W *(RPM/2*60); 
% Calls the KPI function
%KPIs = CalculateKPIs(W, mfr_fuel, LHV, p, mfr_CO2, mfr_NOx);
%=======
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


%>>>>>>> 51126a27191a1963030fa9c8e65cc87d3e1a906e

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

%<<<<<<< HEAD
%% Ideal Diesel Cycle
% Extract initial conditions from actual data
% Find index corresponding to Intake Valve Closure (IVC)
IVC_angle = ValveEvents.CaIVC;  % Crank angle for IVC
P1 = p_filtered_avg(idx_IVC);  % Pressure at IVC
V1 = V_avg(idx_IVC);           % Volume at IVC
% Assume inlet temperature or estimate based on conditions
% <<<<<<< HEAD:ExampleDataSet/Simple.m
% T1 = 295 ;  % K assumed to be ambient
% =======
% T1 = 298.15;  % K (Atmospheric Conditions)
% >>>>>>> 475e35b828e24d257af4f15a6e0772a362edca83:ExampleDataSet/Main.m
%=======
%>>>>>>> 51126a27191a1963030fa9c8e65cc87d3e1a906e


