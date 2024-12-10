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
O2_perc = 14.42; % percentage at exhaust (hardcoded)


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
ValveEvents.CaSOI = -3.2;  %Start of Injection

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

%% Stoichiometric calculations
[stoich_coeffs, reaction_eq, AFR_stoich] = StoichiometricCombustion(fuel_name, SpS, El);

%% Calculate mass flow of air:
avg_O2_load35 = mean(O2_percent_load(4:6));  % Load relevant exhaust data from processed excel
O2_percent_vector = avg_O2_load35 * size(1, 100); %set the correct size for the input
mfr_air = CalculateMassFlowAir(O2_percent_vector, mfr_fuel, AFR_stoich);

%% Calculate Cp and gamma
%[cp, gamma] = calc_cp_gamma(LHV, mfr_fuel, mfr_air);

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
C_p =  1101.6; %Cp manually plugged in from the results of cp and gamma calculations [J/g*K]
gamma = 1.312562;  % gamma manually plugged in from the results of cp and gamma calculations
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


%% Key performance indicators - this still needs to be implemented properly!!
%mfr_CO2 = 
%mfr_Nox = 

% % Power calculation
% Power = W *(RPM/2*60); 
% Calls the KPI function
%KPIs = CalculateKPIs(W, mfr_fuel, LHV, p, mfr_CO2, mfr_NOx);

% Power = W*(RPM/2*60);
% 
% % Calls the KPI function
% KPIs = CalculateKPIs(W, mfr_fuel, LHV, P, mfr_CO2, mfr_NOx);


% %% Key performance indicators
% % % Power calculation
Power = W*(RPM/2*60);
% 
% % Calls the KPI function
% %Example data for diesel
% x = 12;
% LHV = 43e6; %MJ/kg
% mean_mfr_fuel = mean(mean(mfr_fuel,1))%Temporary mfr_fuel values for each CA experiment until problem solved 
% MW_Fuel = 200; %Molar weight of fuel
% % Calls the KPI function
% KPIs = CalculateKPIs(W, LHV, P, mfr_air, x, NOx_ppm_CA, MW_Fuel, mean_mfr_fuel)
% crankAngles = (15:21)'; % Assuming crank angles are 15 to 21 degrees
% fprintf('Crank Angle (°)\tbsNOx\n');
% fprintf('----------------\t-----\n');
% for i = 1:length(KPIs.bsNOx)
%     fprintf('%15.2f\t%5.2f\n', crankAngles(i), KPIs.bsNOx(i));
% end

% %% Given variables
% % Define mole fractions for each species
% % fractions represent typical air composition
% O2_frac = 0.1447;      % Oxygen fraction
% CO2_frac = 0.0467;    % Carbon dioxide fracitons
% N2_frac = 0.7808;      % Nitrogen fraction
% % Calculate water vapor fraciton by subtraction
% H2O_frac = 0.2095 - O2_frac - CO2_frac;
% CO_frac = 0;      % Carbon monoxide fraction
% 
% % Combustion and thermal parameters
% % These should also be read from a table/structure.
% m_dot_fuel = mean(mfr_fuel);      % Mass flow rate of fuel (kg/s)
% Q_dot = Q_combustion_percycle;    % Heat transfer rate (W)
% T_initial = T_int;       % Initial temperature (K)
% tolerance = 1e2;          % Acceptable error in heat transfer (W)
% deltaT = delta_T_percycle;             % Initial guess for temperature change (K)
% error = Inf;              % Initialize error to infinite
% m_dot_air = mfr_air;
% m_dot_tot = m_dot_air + m_dot_fuel; % Total mass flow rate
% 
% %% Find species
% % Locate indexes of specific species in the database
% iSp = myfind({Sp.Name}, {'O2', 'CO2', 'N2','H2O'}); % Find indexes of these species
% SpS = Sp(iSp); % Create subset of species based on found indexes
% NSp = length(SpS); % Number of species
% 
% %% Convert mole fractions to mass fractions
% 
% % Collect mole fractions into an array
% moleFractions = [O2_frac, CO2_frac, N2_frac, H2O_frac];
% 
% % Get molar masses of species
% Mi = [SpS.Mass]; % Molar masses in kg/mol
% 
% % Calculate mass composition based on mole fractions and molar masses
% massComposition = (moleFractions .* Mi) / sum(moleFractions .* Mi);
% 
% % Initialize counter for iteration tracking
% counter = 0;
% cpCum = []; % Array to potentially store cumulative cp values
% TCum = [];  % Array to potentially store cumulative temperatures
% 
% %% Iterative calculation to find temperature change
% % Iterate until calculated heat transfer is within tolerance of given heat transfer
% while error > tolerance
%     % Calculate new temperature
%     T = T_initial + deltaT;
% 
%     % Compute average temperature for cp calculation
%     T_avg = (T + T_initial) / 2;
% 
%     % Calculate specific heat capacity at constant pressure
%     cp = compute_cp(T_avg, SpS, massComposition);
% 
%     % Calculate heat transfer based on current cp and temperature change
%     Q_dot_calculated = m_dot_tot * cp * deltaT;
% 
%     % Compute the error between calculated and given heat transfer
%     error = abs(Q_dot_calculated - Q_dot);
% 
%     % Increment temperature change and iteration counter
%     deltaT = deltaT + 0.1;
%     counter = counter + 1;
% 
%     % for plotting
%     TCum = [TCum T_avg];
%     cpCum = [cpCum cp];
% end
% 
% % Display success message and results
% disp("Great Success, cp = ");
% disp(cp);
% 
% plot(TCum,cpCum)
% xlabel("Temperature")
% ylabel("C_p")
% title("c_p for post-combustion mixture vs temperature")
% subtitle("linear behaviour?")
% 
% %% Repeated display of results
% disp("Great Success, cp = ");
% disp(cp);
% disp("Counter = ")
% disp(counter)
% disp("DeltaT = ")
% disp(deltaT)
% 
% %% Find temperature where specific heat capacity becomes negative
% % Start at a high temperature
% T_test = 5000;
% cp_test = 5;
% 
% % Decrease temperature until cp becomes negative
% while cp_test > 0
%     cp_test = compute_cp(T_test,SpS,massComposition);
%     T_test = T_test + 1;
% end
% 
% % Compute cp at a temperature slightly below the negative cp point
% cp_show = compute_cp(T_test-2000,SpS,massComposition);
% 
% % Display results of cp investigation
% fprintf("c_p for T > %i is negative! = %f \n", T_test, cp_test);
% fprintf("c_p for %i is %f \n", T_test-2000, cp_show);
% 
% %% Compute heat capacity ratio (gamma)
% % Calculate ratio of cp to cv at maximum combustion temperature
% gamma = compute_cp(T,SpS,massComposition) / compute_cv(T,SpS,massComposition);

%% Calculate aROHR
 
aROHR_avg = aROHR(p_filtered_avg, V_avg, resolution, gamma);

% Plot the apparent Rate of Heat Release
figure;
plot(Ca(:, 1), aROHR_avg, 'LineWidth', 1.5);
xlabel('Crank Angle (°)');
ylabel('Apparent Rate of Heat Release [J/°]');
title('Apparent Rate of Heat Release (Average)');
grid on;

%% Calculate Apparent Heat Release
aHR_avg = aHR(aROHR_avg, resolution);  % Assuming aHR function is already defined

% Plot the Apparent Heat Release
figure;
plot(Ca(:, 1), aHR_avg, 'LineWidth', 1.5);  % Ca(:,1) is the crank angle array
xlabel('Crank Angle (°)');
ylabel('Apparent Heat Release [J]');
title('Apparent Heat Release (Average)');
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
IVC_angle = ValveEvents.CaIVC;  % Crank angle for IVC
P1 = p_filtered_avg(idx_IVC);  % Pressure at IVC
V1 = V_avg(idx_IVC);           % Volume at IVC
% Assume inlet temperature or estimate based on conditions
% <<<<<<< HEAD:ExampleDataSet/Simple.m
% T1 = 295 ;  % K assumed to be ambient
% =======
% T1 = 298.15;  % K (Atmospheric Conditions)




