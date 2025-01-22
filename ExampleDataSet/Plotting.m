%% ============================================
% MATLAB Script: Plot Comparison and P-V Diagrams for Most Efficient Experiments
% ============================================
% This script groups experiments by composition, identifies the most efficient 
% experiment for each composition (based on calculated efficiency), and plots 
% its P-V diagram along with aHR and aROHR curves. It dynamically computes 
% the specific heat ratio (γ) using engine and exhaust parameters.
% ============================================

%% Load NASA Data for Species Structure
global Runiv
Runiv = 8.314;  % Universal gas constant [J/(mol·K)]

% Load specified species from NASA thermal database. Diesel will act as a 
% surrogate for various fuel types (HVO, GTL, etc.)
[SpS, El] = myload('Nasa/NasaThermalDatabase.mat', {'N2', 'O2', 'CO2', 'H2O', 'Diesel'});

%% Check if Table T Exists
if ~exist('T', 'var')
    error('Table T does not exist in the workspace. Please run the data integration script first.');
end
disp(T(1:5,:));  % Display first few rows for verification

%% Define Cylinder Parameters
Cyl.Bore = 104;       % mm
Cyl.Stroke = 85;      % mm
Cyl.CompressionRatio = 21.5;
Cyl.ConRod = 136.5;   % mm
Cyl.TDCangle = 180;   % degrees

%% Define Lower Heating Values
LHV = struct();
LHV.Diesel = 43e3;
LHV.HVO = 44.29e3;
LHV.GTL = 44e3;
%% Constants 
C_p = 1101.6;         % [J/g*K]
% We'll calculate gamma dynamically, so no constant gamma here.

%% Define Crank Angle Resolution and Cycle Parameters
resolution = 0.2;  % Degrees crank angle resolution
NdatapointsPerCycle = floor(720 / resolution); % Data points per cycle
crank_angle = (0:NdatapointsPerCycle-1) * resolution'; % 0 to 720 degrees


fprintf('Defined crank angle resolution: %.1f degrees.\n', resolution);
fprintf('Number of data points per cycle: %d.\n', NdatapointsPerCycle);

%% Predefine Engine Parameters for γ Calculation
% These may be adjusted or determined dynamically per experiment
AFR = 14.7;           % Stoichiometric air-fuel ratio for gasoline (adjust as needed)
defaultRPM = 1500;    % Default RPM if not available in metadata
defaultMFR = 0.15;
%% Initialize Arrays for Efficiency and Composition Grouping
numGroups = height(T);
efficiency_all = zeros(numGroups,1);
composition_all = cell(numGroups,1);
aHR_all = cell(numGroups,1);
aROHR_all = cell(numGroups,1);

%% Loop Through All Groups in Table T to Compute Metrics
for rowIdx = 1:numGroups
    % Extract group information
    uniqueID = T.UniqueID{rowIdx};
    metadata = T.Metadata(rowIdx);  % Assuming metadata is a struct in each row
    averageData = T.AverageCycleData{rowIdx};
    
    fprintf('Processing Group %d/%d: %s\n', rowIdx, numGroups, uniqueID);
    
    % Extract averaged variables
    avgPressure = averageData.AvgPressure;  
    avgMassFlow = averageData.AvgMassFlow;  
    avgCurrent = averageData.AvgCurrent;    
    
    mfr_fuel = CalculateMassFlowFuel(avgMassFlow, avgCurrent, crank_angle', defaultRPM, 0.5);
    if mfr_fuel < 0.07
        warning('Error in Mass flow sensor data using default value for mfr')
        mfr_fuel = defaultMFR;
    end
    % Calculate cylinder volume for each crank angle [m^3]
    volume = CylinderVolume(crank_angle, Cyl) * 1e-9; % convert mm^3 to m^3
    
    % Retrieve additional data for exhaust composition
    additionalData = T.AdditionalData{rowIdx};
    if ~isempty(additionalData) && isstruct(additionalData)
        O2percent = additionalData.O2;
        CO2percent = additionalData.CO2;
    else
        % Default values if additional data is missing or invalid
        O2percent = 0;
        CO2percent = 0;
    end

    % Use RPM from metadata if available, otherwise fallback to defaultRPM
    if isfield(metadata, 'RPM')
        RPM = metadata.RPM;
    else
        RPM = defaultRPM;
    end
    
    C = metadata.C;
    FuelType = metadata.FuelType;
    % Compute dynamic gamma over the crank angle using CalculateGamma
    try
        gamma_curve = CalculateGamma(SpS, volume', avgPressure, O2percent, CO2percent, mfr_fuel, AFR, RPM);
    catch ME
        warning('Failed to calculate gamma for group %s: %s', uniqueID, ME.message);
        gamma_curve = ones(size(avgPressure)) * 1.312562;  % fallback value
    end
    
    % Compute aROHR and aHR for the current experiment using dynamic gamma
    aROHR_result = get_aROHR(avgPressure, volume', gamma_curve); 
    aHR_result = get_aHR(aROHR_result); 
 
    % Store the results for later use
    aHR_all{rowIdx} = aHR_result; 
    aROHR_all{rowIdx} = aROHR_result; 

    % Verify data lengths
    if length(avgPressure) ~= length(crank_angle) || ...
       length(avgMassFlow) ~= length(crank_angle) || ...
       length(avgCurrent) ~= length(crank_angle)
        warning('Mismatch in data lengths for group %s. Skipping plotting.', uniqueID);
        continue;
    end
    % compute the composition of LHV in the current fuel mix
    CompLHV= C/100*LHV.(FuelType) + (1-C/100)*LHV.Diesel;
    % Compute net work over the cycle: ∑ P * dV
    work_total = trapz(volume, avgPressure);
    
    Power = work_total * (RPM / (2 * 60)) ;
    
    
    
    % Calculate efficiency as the ratio of net work to total heat release
    efficiency_all(rowIdx) = Power / (mfr_fuel * CompLHV)
    
    % Create a composition key based on metadata fields
    composition_all{rowIdx} = sprintf('C%.2f_Fuel%s', metadata.C, metadata.FuelType);
end

%% Group by Composition and Select Most Efficient Experiment per Composition
unique_compositions = unique(composition_all);
most_efficient_idx = zeros(length(unique_compositions),1);

for k = 1:length(unique_compositions)
    comp = unique_compositions{k};
    indices = find(strcmp(composition_all, comp));
    [~, maxLocalIdx] = max(efficiency_all(indices));
    most_efficient_idx(k) = indices(maxLocalIdx);
end

selected_indices = most_efficient_idx;
colors = lines(length(selected_indices));

%% Plot aHR for Most Efficient Experiments
figure;
hold on;
for i = 1:length(selected_indices)
    idx = selected_indices(i);
    uniqueID = T.UniqueID{idx};
    % Subtract 360 from crank_angle for plotting
    plot(crank_angle - 360, aHR_all{idx}, 'LineWidth', 1.5, 'Color', colors(i, :), ...
         'DisplayName', sprintf('%s', uniqueID));
end
xlabel('Crank Angle (°)');
ylabel('Apparent Heat Release [J]');
title('aHR for Most Efficient Experiments by Composition');
legend('show');
grid on;
xlim([-10,50]);
hold off;

%% Plot aROHR for Most Efficient Experiments
figure;
hold on;
for i = 1:length(selected_indices)
    idx = selected_indices(i);
    uniqueID = T.UniqueID{idx};
    % Subtract 360 from crank_angle for plotting
    plot(crank_angle - 360, aROHR_all{idx}, 'LineWidth', 1.5, 'Color', colors(i, :), ...
         'DisplayName', sprintf('%s', uniqueID));
end
xlabel('Crank Angle (°)');
ylabel('Apparent Rate of Heat Release [J/°]');
title('aROHR for Most Efficient Experiments by Composition');
legend('show');
grid on;
xlim([-10,50]);
hold off;

%% Plot P-V Diagram for Each Most Efficient Experiment
for j = 1:length(selected_indices)
    rowIdx = selected_indices(j);
    uniqueID = T.UniqueID{rowIdx};
    averageData = T.AverageCycleData{rowIdx};
    
    % Re-extract data and recompute volume for the selected experiment
    avgPressure = averageData.AvgPressure;
    crank_angle = (0:NdatapointsPerCycle-1) * resolution;
    volume = CylinderVolume(crank_angle, Cyl) * 1e-9; % convert mm^3 to m^3
    
    figure('Name', sprintf('P-V Diagram - Group %s', uniqueID), ...
           'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
    plot(volume, avgPressure, 'm', 'LineWidth', 1.5);
    title(sprintf('Pressure-Volume (P-V) Diagram\nGroup: %s', uniqueID));
    xlabel('Volume (m^3)');
    ylabel('Pressure (Pa)');
    grid on;
    legend('P-V Curve', 'Location', 'best');
end

fprintf('Completed plotting for most efficient experiments of each composition.\n');

%% ===========================================
% End of Script
% ===========================================
