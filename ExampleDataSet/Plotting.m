% ============================================
% MATLAB Script: Plot Comparison and P-V Diagrams for Most Efficient Experiments
% ============================================
% This script groups experiments by composition, identifies the most efficient 
% experiment for each composition (based on calculated efficiency), and plots 
% its P-V diagram along with aHR and aROHR curves.
% ============================================

%% Check if Table T Exists
if ~exist('T', 'var')
    error('Table T does not exist in the workspace. Please run the data integration script first.');
end
T
%% Define Cylinder Parameters
Cyl.Bore = 104; % mm
Cyl.Stroke = 85; % mm
Cyl.CompressionRatio = 21.5;
Cyl.ConRod = 136.5; % mm
Cyl.TDCangle = 180; % degrees

%% Constants 
C_p =  1101.6; % [J/g*K]
gamma = 1.312562;  % 


%% Define Crank Angle Resolution and Cycle Parameters
resolution = 0.2;  % Degrees crank angle resolution
NdatapointsPerCycle = floor(720 / resolution); % Data points per cycle
crank_angle = (0:NdatapointsPerCycle-1) * resolution; % 0 to 720 degrees
crank_angle_trimmed = crank_angle(1:end);

fprintf('Defined crank angle resolution: %.1f degrees.\n', resolution);
fprintf('Number of data points per cycle: %d.\n', NdatapointsPerCycle);

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
    
    % Calculate cylinder volume for each crank angle
    volume = CylinderVolume(crank_angle, Cyl); % cubic mm
    
    T.AdditionalData{rowIdx}
    % Compute aROHR and aHR for the current experiment
    aROHR_result = get_aROHR(avgPressure, volume', gamma); 
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
    
    % Compute net work over the cycle: ∑ P * dV
    dV = diff(volume);
    work_total = sum(avgPressure(1:end-1) .* dV(:));
    
    % Total heat released is the final value of aHR
    Q_total = aHR_result(end);
    
    % Calculate efficiency as the ratio of net work to total heat release
    efficiency_all(rowIdx) = work_total / Q_total;
    
    % Create a composition key based on metadata fields
    composition_all{rowIdx} = sprintf('C%.2f_Fuel%s', metadata.C, metadata.FuelType);
end

%% Group by Composition and Select Most Efficient Experiment per Composition
unique_compositions = unique(composition_all)
most_efficient_idx = zeros(length(unique_compositions),1);

for k = 1:length(unique_compositions)
    comp = unique_compositions{k};
    indices = find(strcmp(composition_all, comp));
    [~, maxLocalIdx] = max(efficiency_all(indices));
    most_efficient_idx(k) = indices(maxLocalIdx);
end

selected_indices = most_efficient_idx
colors = lines(length(selected_indices));

%% Plot aHR for Most Efficient Experiments
figure;
hold on;
for i = 1:length(selected_indices)
    idx = selected_indices(i);
    uniqueID = T.UniqueID{idx};
    plot(crank_angle_trimmed, aHR_all{idx}, 'LineWidth', 1.5, 'Color', colors(i, :), ...
         'DisplayName', sprintf('%s', uniqueID));
end
xlabel('Crank Angle (°)');
ylabel('Apparent Heat Release [J]');
title('aHR for Most Efficient Experiments by Composition');
legend('show');
grid on;
xlim([-45, 135]);
hold off;

%% Plot aROHR for Most Efficient Experiments
figure;
hold on;
for i = 1:length(selected_indices)
    idx = selected_indices(i);
    uniqueID = T.UniqueID{idx};
    plot(crank_angle_trimmed, aROHR_all{idx}, 'LineWidth', 1.5, 'Color', colors(i, :), ...
         'DisplayName', sprintf('%s', uniqueID));
end
xlabel('Crank Angle (°)');
ylabel('Apparent Rate of Heat Release [J/°]');
title('aROHR for Most Efficient Experiments by Composition');
legend('show');
grid on;
xlim([-45, 135]);
hold off;

%% Plot P-V Diagram for Each Most Efficient Experiment
for j = 1:length(selected_indices)
    rowIdx = selected_indices(j);
    uniqueID = T.UniqueID{rowIdx};
    averageData = T.AverageCycleData{rowIdx};
    
    % Re-extract data and recompute volume for the selected experiment
    avgPressure = averageData.AvgPressure;
    crank_angle = (0:NdatapointsPerCycle-1) * resolution;
    volume = CylinderVolume(crank_angle, Cyl);
    
    figure('Name', sprintf('P-V Diagram - Group %s', uniqueID), ...
           'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
    plot(volume, avgPressure, 'm', 'LineWidth', 1.5);
    title(sprintf('Pressure-Volume (P-V) Diagram\nGroup: %s', uniqueID));
    xlabel('Volume (cubic mm)');
    ylabel('Pressure (kPa)');
    grid on;
    legend('P-V Curve', 'Location', 'best');
end

fprintf('Completed plotting for most efficient experiments of each composition.\n');
