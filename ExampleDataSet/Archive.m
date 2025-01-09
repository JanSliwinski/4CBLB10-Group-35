%% Archive
% Contains parts archived from the main code due to them becoming
% unnecesarily at some point
%

%% Making Path for files of different blends ( archived on: 31/12/2024 by Kata)

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

%% Valve Events (Crank Angles)     (archived on: 31/12/2024 by Kata, reason: important valve events defined in another way)
ValveEvents.CaIVO = -355;
ValveEvents.CaIVC = -135;
ValveEvents.CaEVO = 149;
ValveEvents.CaEVC = -344;
ValveEvents.CaSOI = -3.2;  % Start of Injection

%% Load and Reshape Data            (archived on: 31/12/2024 by Kata, reason: changed data loading to use table T)
dataFileName = fullfile('Data' , 'processed_Data_experiment1_load3.5.txt');
dataIn = table2array(readtable(dataFileName));

%% Filter Pressure Data             (archived on: 1/1/2025 by Kata, reason: table T already contains filtered pressure data so this became unnecesary)
polynomialOrder = 3;
frameLength = 21;  % Must be odd

% Initialize the filtered pressure matrix
p_filtered = zeros(size(p));

% Apply the filter to each column
for i = 1:Ncycles
    p_filtered(:, i) = SGFilter(p(:, i), polynomialOrder, frameLength, 0);
end

disp('Data filtered and reshaped into cycles');

%% load the excelfile               (archived on: 1/1/2024 by Kata, reason: exhaust information is now loaded in from table T)
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


%% Code from GenerateKPITable.m
 % Reshape Data
        resolution = 0.2;  % Degrees crank angle resolution
        n_datapoints_per_cycle = 720 / resolution;
        n_cycles = size(data_in, 1) / n_datapoints_per_cycle;

        ca = reshape(data_in(:, 1), [], n_cycles);          % Crank angle in degrees
        p = reshape(data_in(:, 2), [], n_cycles) * 1e5;     % Pressure in Pa
        mfr_fuel = true_mfr_fuel;    % Fuel mass flow rate (kg/s)

        % Apply Savitzky-Golay filter to pressure data
        polynomial_order = 3; % Adjust based on noise level
        frame_length = 21;    % Must be odd
        
        % Initialize the filtered pressure matrix
        p_filtered = zeros(size(p));
        
        % Apply the filter to each cycle
        for j = 1:n_cycles
            p_filtered(:, j) = SGFilter(p(:, j), polynomial_order, frame_length, 0);
        end
        
        % Replace raw pressure data with filtered data for further calculations
        p = p_filtered;
 %% Calculate mass flow of air:          (archived on 05/1/2025 by Kata, reason: was already commented out, plus air flow rate is calculated elsewhere now)

% mfr_air = CalculateMassFlowAir(O2_percent, mfr_fuel, AFR_stoich);
% mfr_air = CalculateMassFlowAir(O2_percent, true_mfr_fuel, AFR_stoich);
% fprintf('Mass flow rate for air: %.6f g/s\n', mfr_air)

%% Calculate Heat of combustion and Temperature at exhaust - LHV way                (archived on 5/1/2025 by Kata, reason: these calculations ar enot used anymore)
% [T_exh, Q_combustion_LHV, m_combusted] = Calc_Q_LHV(C_p, mfr_fuel, mfr_air, RPM, W, LHV, T_int);

%% Calculate Heat of combustion and Temperature at exhaust - aROHR way              (archived on 5/1/2025 by Kata, reason: these calculations ar enot used anymore)
%[Q_combustion_aROHR] = Calc_Q_aROHR(p_filt, V_all, Ca, ValveEvents, gamma);
% disp(['Q combustion aROHR way: ', num2str(Q_combustion_aROHR), ' J']);

%% Thermal efficiency of the engine                                                 (archived on 5/1/2025 by Kata, reason: these calculations ar enot used anymore)
% efficiency_LHV = (W / Q_combustion_LHV) *100; % efficiency for each cycle
% disp(['Calculated average thermal efficiency(LHV): ', num2str(efficiency_LHV), ' %']);
% 
% efficiency_LHV = (W / Q_combustion_aROHR) *100; % efficiency for each cycle
% disp(['Calculated average thermal efficiency(aROHR): ', num2str(efficiency_LHV), ' %']);

%% Calculate thermodynamic properties for each cycle - THIS STILL NEEDS TO BE IMPLEMENTED PROPERLY - NEED TO CALCULAT T_EXHAUST FOR IT SOMEHOW 
intake_species = [2, 3];           % Example species (O2 and N2)                    (archived on 5/1/2025 by Kata, reason: was never implemented)
exhaust_species = [4, 5, 3];       % Example species (CO2, H2O, and N2)
Y_int = [0.21, 0.79];              % Mole fractions for intake
Y_exh = [0.12, 0.18, 0.70];        % Mole fractions for exhaust

% % Call the function
%[Delta_H_all, Delta_U_all, Delta_S_all] = ThermoProperties(T_int, T_exh, SpS, Ncycles, Ca, intake_species, exhaust_species, Y_int, Y_exh);
 
% % Calculate averages
% Delta_H_avg = mean(Delta_H_all, 2);
% Delta_U_avg = mean(Delta_U_all, 2);
% Delta_S_avg = mean(Delta_S_all, 2);

%% Calculate mass flow of fuel; Bart's method
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