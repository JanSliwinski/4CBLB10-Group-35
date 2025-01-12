% %% Initialization
% clear; clc; close all;
% 
% % Load your Table T (from "T.mat")
% data = load("T.mat");
% T = data.T;

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
M_HVO = 226;    % Molar mass (g/mol)
Density_HVO = 0.78;
Density_Diesel = 0.85;
LHV_diesel = 43e3;  %Lower heating value given in the project guide for Diesel B7 J/g
Ca = 0:0.2:720-0.2;
% LHV_HVO = ?

%% Load NASA Data (if needed)
global Runiv
Runiv = 8.314; 

[SpS, El] = myload('Nasa/NasaThermalDatabase.mat', {'N2', 'O2', 'CO2', 'H2O', 'Diesel'});

%% Volume (Engine Geometry)
Cyl.Bore = 104 * mm;               
Cyl.Stroke = 85 * mm;              
Cyl.CompressionRatio = 21.5;       
Cyl.ConRod = 136.5 * mm;           
Cyl.TDCangle = 180;                
volume = CylinderVolume(Ca, Cyl)';
disp('Cylinder volume calculated / cycle');

%% Stoichiometric calculations
[stoich_coeffs, reaction_eq, AFR_stoich] = StoichiometricCombustion('Diesel', SpS, El);
stoich_coeffs.fuel

%% Loop Through All Groups in Table T
numGroups = height(T);
fprintf('Starting to plot data for %d groups.\n', numGroups);

% --- PREINITIALIZE EVERYTHING AS CELL ARRAYS --- %
ID           = cell(numGroups,1);  % Stores uniqueID
aHR_all      = cell(numGroups,1);
aROHR_all    = cell(numGroups,1);
Ca50Idx      = cell(numGroups,1);
Ca50aROHR    = cell(numGroups,1);
Ca50aHR      = cell(numGroups,1);

for rowIdx = 1:numGroups
    try
        % Extract group information
        uniqueID    = T.UniqueID{rowIdx};
        metadata    = T.Metadata(rowIdx, :);
        averageData = T.AverageCycleData{rowIdx};
        
        % Only process if FuelType is 'HVO' & C == 0
        if metadata.FuelType == 'HVO' & metadata.C == 0
            fprintf('Processing Group %d/%d: %s\n', rowIdx, numGroups, uniqueID);
            
            [dataIn, ExhaustData, Ca, p_avg, S_current, mfr_fuel, ...
             CO_percent_load, HC_ppm_load, NOx_ppm_load, ...
             CO2_percent_load, O2_percent_load, lambda_load] = loadingfromT(T, uniqueID, bara);
            p_avg = averageData.AvgPressure;
            
            true_mfr_fuel = mean(mfr_fuel);
    
            gamma = CalculateGamma(SpS, volume, p_avg, O2_percent_load, ...
                                   CO2_percent_load, true_mfr_fuel, AFR_stoich, RPM);
            aROHR = get_aROHR(p_avg, volume, gamma);
            
            % Zero out or isolate region around 355 -> 410 deg CA
            idxStart = 355 / 0.2;
            idxEnd   = (50 + 360) / 0.2;
            aROHR(1:idxStart)    = 0;
            aROHR(idxEnd:end)    = 0;
            aHR = get_aHR(aROHR);
            
            fprintf('Successfully processed group %d: %s\n', rowIdx, uniqueID);

            % --- STORE RESULTS IN CELL ARRAYS ---
            ID{rowIdx}        = uniqueID;   % <--- curly braces
            aHR_all{rowIdx}   = aHR; 
            aROHR_all{rowIdx} = aROHR; 
            
            % CA50
            halfMax = max(aHR)/2; 
            [~, Ca50idx] = min(abs(aHR - halfMax));
            Ca50Idx{rowIdx}   = Ca50idx/5 - 360;  % dividing by 5 because 1 deg = 5 steps if 0.2 deg/step
            Ca50aROHR{rowIdx} = aROHR(Ca50idx);
            Ca50aHR{rowIdx}   = aHR(Ca50idx);
        end
        
    catch ME
        % Log the error with detailed info
        fprintf('\nERROR in group %d (%s):\n', rowIdx, uniqueID);
        fprintf('Error Message: %s\n', ME.message);
        fprintf('Error Location: %s (Line %d)\n', ME.stack(1).name, ME.stack(1).line);
        fprintf('Skipping to next group...\n\n');
        continue;
    end
end

% --- REMOVE EMPTY CELLS & CONVERT TO STRINGS ---
nonEmptyCells = ID(~cellfun(@isempty, ID));
ID = string(nonEmptyCells);

%% Loop to plot aROHR and aHR

a=42;
b=43;


figure;

subplot(1,2,1)
hold on;
% for i = 1:length(aROHR_all)
for i = [a,b]
    if ~isempty(aROHR_all{i})
        % Plot the aROHR curve
        plot(Ca, aROHR_all{i}, 'LineWidth', 0.5, ...
             'DisplayName', sprintf('%s', T.UniqueID{i}));
    end
end
xlabel('Crank Angle (°)');
ylabel('Apparent Rate of Heat Release [J/deg]');
title('aROHR');
legend('show', 'Location','northeast');  % Let MATLAB auto-build the legend
grid on;
xlim([-10,50]);
hold off;

subplot(1,2,2)
hold on;
% for i = 1:length(aHR_all)
for i = [a,b]
    if ~isempty(aHR_all{i})
        plot(Ca, aHR_all{i}, 'LineWidth', 0.5, ...
             'DisplayName', sprintf('%s', T.UniqueID{i}));
        % Plot the CA50 marker
        if i == b
            scatter(Ca50Idx{i}, Ca50aHR{i}, 50, 'filled', ...
                    'MarkerFaceColor', 'r', 'DisplayName', 'CA50');
        else
            scatter(Ca50Idx{i}, Ca50aHR{i}, 50, 'filled', ...
                    'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
        end
    end
end
xlabel('Crank Angle (°)');
ylabel('Apparent Heat Release [J]');
title('aHR');
legend('show', 'Location','southeast');
grid on;
xlim([-10,50]);
hold off;

sgtitle("Diesel Injection Angle Optimization")

fprintf('Completed plotting for all groups.\n');
