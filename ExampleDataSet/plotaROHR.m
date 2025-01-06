%% Initialization
clear all
if ~exist('T', 'var')
    try
        load('T.mat', 'T');
    catch ME
        fprintf('Run Data Loading file!\n');
        % Consider adding error handling here
        rethrow(ME);
    end
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
Ca = 0:0.2:720-0.2;
% LHV_HVO = ?

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
volume = CylinderVolume(Ca,Cyl)';
disp('Cylider volume calculated / cycle');

%% Stoichiometric calculations
[stoich_coeffs, reaction_eq, AFR_stoich] = StoichiometricCombustion('Diesel', SpS, El);
stoich_coeffs.fuel

%% Loop Through All Groups in Table T
numGroups = height(T);
fprintf('Starting to plot data for %d groups.\n', numGroups);
ID = [];

for rowIdx = 1:numGroups
    try
        % Extract group information
        uniqueID = T.UniqueID{rowIdx};
        metadata = T.Metadata(rowIdx, :);
        averageData = T.AverageCycleData{rowIdx};
        
        fprintf('Processing Group %d/%d: %s\n', rowIdx, numGroups, uniqueID);
        
        [dataIn, ExhaustData, Ca, p_filt, S_current, mfr_fuel, ...
         CO_percent_load, HC_ppm_load, NOx_ppm_load, ...
         CO2_percent_load, O2_percent_load, lambda_load] = loadingfromT(T, uniqueID, bara);
        
        true_mfr_fuel = mean(mfr_fuel);
        p_filtRough = sgolayfilt(p_filt, 2,15);

        gamma = CalculateGamma(SpS,volume,p_filtRough,O2_percent_load,CO2_percent_load,true_mfr_fuel,AFR_stoich,RPM);
        % Calculate ROHR and HR
        aROHR = get_aROHR(p_filtRough, volume, gamma);
        
        % Isolate peak
        idxStart = 355 / 0.2;
        idxEnd = (50 + 360) / 0.2;
        aROHR(1:idxStart) = 0;
        aROHR(idxEnd:end) = 0;
        aHR = get_aHR(aROHR);
        
        fprintf('Successfully processed group %d: %s\n', rowIdx, uniqueID);
        ID{rowIdx} = uniqueID;
        
    catch ME
        % Log the error with detailed information
        fprintf('\nERROR in group %d (%s):\n', rowIdx, uniqueID);
        fprintf('Error Message: %s\n', ME.message);
        fprintf('Error Location: %s (Line %d)\n', ME.stack(1).name, ME.stack(1).line);
        fprintf('Skipping to next group...\n\n');
        continue;
    end

    % Store the results
    aHR_all{rowIdx} = aHR; % Store aHR for each file
    aROHR_all{rowIdx} = aROHR; %Store aROHR for each file
end

nonEmptyCells = ID(~cellfun(@isempty, ID));
ID = string(nonEmptyCells);

%% Loop to plot aROHR
figure;
hold on;
for i = 1:length(aROHR_all)
    if ~isempty(aROHR_all{i})
        uniqueID = T.UniqueID{i};
        plot(Ca, aROHR_all{i}, 'LineWidth', 1.5, 'DisplayName', sprintf(' %d', uniqueID));
    end
end
xlabel('Crank Angle (°)');
ylabel('Apparent Heat Release [J]');
title('aROHR for All Files');
legend(ID');
grid on;
xlim([-10,50]);
hold off;
fprintf('Completed plotting for all groups.\n');

%% Loop to plot aHR
figure;
hold on;
for i = 1:length(aHR_all)    % Changed from aROHR_all to aHR_all
    if ~isempty(aHR_all{i})  % Changed from aROHR_all to aHR_all
        uniqueID = T.UniqueID{i};
        plot(Ca, aHR_all{i}, 'LineWidth', 1.5, 'DisplayName', sprintf(' %d', uniqueID));  % Changed from aROHR_all to aHR_all
    end
end
xlabel('Crank Angle (°)');
ylabel('Cumulative Heat Release [J]');  % Changed ylabel to reflect aHR instead of aROHR
title('aHR for All Files');            % Changed title to reflect aHR
legend(ID');
grid on;
xlim([-10,50]);
hold off;
fprintf('Completed plotting for all groups.\n');