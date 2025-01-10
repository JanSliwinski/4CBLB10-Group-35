warning off
%% Initialization
if ~exist("T","var")    % make sure that table T is defined
    load('T.mat','T')
end
clearvars -except T; close all; clc % clean out the workspace and the command window
addpath("Functions", "Nasa");savepath  % Add necessary paths for Nasa tables

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

ID = 'L50I20C100FuelGTL'; %DEFINE ID OF THE EXPERIMENT DATA YOU WANT TO LOAD IN!
%run fucntion to load in all relevant data
[dataIn, ExhaustData, Ca, p_avg, S_current, mfr_fuel, CO_percent_load, HC_ppm_load, NOx_ppm_load, CO2_percent_load, O2_percent_load, lambda_load] = loadingfromT(T, ID, bara);
IDsforKPI = table({'L50I14C100FuelGTL', 'L50I16C100FuelGTL', 'L50I18C100FuelGTL', 'L50I20C100FuelGTL'}');

%% Define Fuel used and applicable LHV - CHANGE THE LINES IN THIS SECTION IF RUN WITH A DIFFERENT FUEL!!!
fuel_used = 'GTL100';
perc_blend = 1; %fraction of the blended in fuel (HVO or GTL) should be set to 0 if you are using diesel
x_blend = x_diesel;  %carbon atoms in the given fuel, can be: x_diesel, x_HVO or x_GTL
LHV_blend = LHV_GTL; %LHV of the given fuel, can be: LHV_diesel, LHV_HVO or LHV_GTL


%% Calculate LHV and x for the fuel used
perc_diesel = 1-perc_blend;     % percentage of diesel in the fuel used 
LHV = LHV_blend * perc_blend + LHV_diesel * perc_diesel;    % calculate the LHV for the fuel
x = x_blend * perc_blend + x_diesel * perc_diesel;  % calculate the number of carbon atoms in the fuel
MW_fuel = M_diesel;     % Molar mass of fuel - keep this on M_diesel

%% Load NASA Data (if needed)
global Runiv
Runiv = 8.314;
[SpS, El] = myload('Nasa/NasaThermalDatabase.mat', {'N2', 'O2', 'CO2', 'H2O', 'Diesel'}); 

%% Compute massflow
sps_fuel_name = 'Diesel';   % fuel name for Sps calculations (keep it on diesel!)
if mean(mfr_fuel) > 0.1 || mean(mfr_fuel) < 0.3     % check if measured mfr fuel is realistic
    mfr_fuel = mean(mfr_fuel);  % set mfr fuel to the mean of all measured data
    [~,~,AFR_stoich] = StoichiometricCombustion(sps_fuel_name, SpS, El); % perform stoichiometric calculation
    mfr_air = CalculateMassFlowAir(O2_percent_load,mfr_fuel,AFR_stoich); % calculate mfr of air
    AFR = mfr_air / mfr_fuel;   % calculate AFR
else
    % Calculate the mass flow rate of fuel using Bart's method.
    % This method is used if the measured mass flow rate of fuel is deemed unrealistic.

    % Constants
    M_C = 12;           % Molar mass of Carbon (g/mol)
    M_H = 1;            % Molar mass of Hydrogen (g/mol)
    M_CO2 = 44;         % Molar mass of CO2 (g/mol)
    mass_CH_ratio = 2;  % Typical Carbon-to-Hydrogen mass ratio for GTL (Gas-to-Liquid fuel)
    CO2_mass_flow_rate = 0.5; % Assumed mass flow rate of CO2 (g/s)

    % Calculate the molar Carbon-to-Hydrogen ratio
    molar_CH_ratio = mass_CH_ratio * (M_H / M_C); % Convert mass CH ratio to molar CH ratio

    % Determine the Hydrogen-to-Carbon ratio (y) for the given fuel
    y = molar_CH_ratio * x; % x represents the Carbon content, given earlier in the code

    % Calculate the molar mass of the fuel (CxHy) in g/mol
    fuel_molar_mass = x * M_C + y * M_H;

    % Convert the mass flow rate of CO2 to a molar flow rate (mol/s)
    CO2_molar_flow_rate = CO2_mass_flow_rate / M_CO2;

    % Determine the molar flow rate of the fuel (CxHy)
    fuel_molar_flow_rate = CO2_molar_flow_rate / x;

    % Calculate the mass flow rate of the fuel (g/s)
    mfr_fuel = mean(fuel_molar_flow_rate * fuel_molar_mass);
end

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

%% Plot Pressure vs. Crank Angle for All Cycles
figure;
plot(Ca, p_avg / bara, 'LineWidth', 1);
xlabel('Crank Angle (°)');
ylabel('Pressure (bar)');
xlim([-360, 360]);
ylim([0, 50]);
title('Pressure vs. Crank Angle for All Cycles');
grid on;

% Highlight a specific cycle
iselect = 10;
hold on;
plot(Ca(iselect), p_avg(iselect) / bara, 'r', 'LineWidth', 2);

% Plot valve events
% YLIM = ylim;
% line([ValveEvents.CaIVC, ValveEvents.CaIVC], YLIM, 'Color', 'b', 'LineWidth', 1);
% line([ValveEvents.CaEVO, ValveEvents.CaEVO], YLIM, 'Color', 'r', 'LineWidth', 1);

set(gca, 'XTick', -360:60:360);
grid on;

%% Calculate Average Volume and Pressure
% V_avg = mean(volume, 2);         % Average volume across all cycles for every CA
% p_avg_num = mean(p_avg, 2);             % Average pressure across all cycles for every CA






%% Calculate Work
W = trapz(volume, p_avg*1e5); % Calculate the area under the averaged p-V curve
disp(['Calculated work: ', num2str(W), ' J']);

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
if exist('O2_percent_load','var')
    gamma = CalculateGamma(SpS,volume,p_avg,O2_percent_load,CO2_percent_load,mfr_fuel,AFR_stoich,RPM);
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

%% Key performance indicators

% Generate KPI table
KPITable = GenerateKPITable(IDsforKPI, mfr_fuel, T, LHV, RPM, AFR_stoich, x, MW_fuel,Cyl, fuel_used);
disp(KPITable)

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
