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

