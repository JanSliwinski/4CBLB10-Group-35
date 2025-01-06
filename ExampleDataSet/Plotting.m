% ============================================
% MATLAB Script: Plot Comparison and P-V Diagrams for All Groups
% ============================================
% This script iterates through all groups in the integrated data table T,
% plotting comparison graphs (Average Pressure, Mass Flow, Current vs Crank Angle)
% and Pressure-Volume (P-V) diagrams using calculated cylinder volumes.
% ============================================

%% Check if Table T Exists
%load("T.mat")
if ~exist('T', 'var')
    error('Table T does not exist in the workspace. Please run the data integration script first.');
end

%% Define Cylinder Parameters
% Units: millimeters (mm) and degrees
Cyl.Bore = 104; % mm
Cyl.Stroke = 85; % mm
Cyl.CompressionRatio = 21.5;
Cyl.ConRod = 136.5; % mm
Cyl.TDCangle = 180; % degrees

%% Constants 

C_p =  1101.6; %Cp manually plugged in from the results of cp and gamma calculations [J/g*K]
gamma = 1.312562;  % gamma manually plugged in from the results of cp and gamma calculations


%% Define Crank Angle Resolution and Cycle Parameters
resolution = 0.2;  % Degrees crank angle resolution
NdatapointsPerCycle = floor(720 / resolution); % Number of data points per cycle

fprintf('Defined crank angle resolution: %.1f degrees.\n', resolution);
fprintf('Number of data points per cycle: %d.\n', NdatapointsPerCycle);

%% Loop Through All Groups in Table T
numGroups = height(T);
fprintf('Starting to plot data for %d groups.\n', numGroups);

for rowIdx = 1:numGroups
    % Extract group information
    uniqueID = T.UniqueID{rowIdx};
    metadata = T.Metadata(rowIdx, :);
    averageData = T.AverageCycleData{rowIdx};
    
    fprintf('Processing Group %d/%d: %s\n', rowIdx, numGroups, uniqueID);
    
    % Extract averaged variables
    avgPressure = averageData.AvgPressure;   % Vector of average pressure values
    avgMassFlow = averageData.AvgMassFlow;   % Vector of average mass flow values
    avgCurrent = averageData.AvgCurrent;     % Vector of average current values
    
    % Define crank angle vector
    crank_angle = (0:NdatapointsPerCycle-1) * resolution; % 0 to 720 degrees
    
    % Calculate cylinder volume for each crank angle
    volume = CylinderVolume(crank_angle, Cyl); % cubic mm
    
    % Verify data lengths
    if length(avgPressure) ~= length(crank_angle) || ...
       length(avgMassFlow) ~= length(crank_angle) || ...
       length(avgCurrent) ~= length(crank_angle)
        warning('Mismatch in data lengths for group %s. Skipping plotting.', uniqueID);
        continue;
    end

    % %% Create Figure for Comparison Graphs and P-V Diagram
    % figure('Name', sprintf('Group %s - Comparison and P-V Diagram', uniqueID), ...
    %        'NumberTitle', 'off', 'Position', [100, 100, 1600, 900]);
    % 
    % % ---- Subplot 1: Average Pressure vs Crank Angle ----
    % subplot(2,2,1);
    % plot(crank_angle, avgPressure, 'b', 'LineWidth', 1.5);
    % title(sprintf('Average Filtered Pressure vs Crank Angle\nGroup: %s', uniqueID));
    % xlabel('Crank Angle (Degrees)');
    % ylabel('Pressure (kPa)'); % Replace 'kPa' with actual units
    % grid on;
    % legend('Avg Pressure', 'Location', 'best');
    % 
    % % ---- Subplot 2: Average Mass Flow vs Crank Angle ----
    % subplot(2,2,2);
    % plot(crank_angle, avgMassFlow, 'r', 'LineWidth', 1.5);
    % title('Average Mass Flow vs Crank Angle');
    % xlabel('Crank Angle (Degrees)');
    % ylabel('Mass Flow (kg/s)'); % Replace 'kg/s' with actual units
    % grid on;
    % legend('Avg Mass Flow', 'Location', 'best');
    % 
    % % ---- Subplot 3: Average Current vs Crank Angle ----
    % subplot(2,2,3);
    % plot(crank_angle, avgCurrent, 'g', 'LineWidth', 1.5);
    % title('Average Current vs Crank Angle');
    % xlabel('Crank Angle (Degrees)');
    % ylabel('Current (A)'); 
    % grid on;
    % legend('Avg Current', 'Location', 'best');
    % 
    % % ---- Subplot 4: Pressure-Volume (P-V) Diagram ----
    % subplot(2,2,4);
    % plot(volume, avgPressure, 'm', 'LineWidth', 1.5);
    % title('Pressure-Volume (P-V) Diagram');
    % xlabel('Volume (cubic mm)');
    % ylabel('Pressure (kPa)'); 
    % grid on;
    % legend('P-V Curve', 'Location', 'best');
    % 
    % % Adjust layout for better visibility
    % sgtitle(sprintf('Comparison and P-V Diagram for Group %s', uniqueID), 'FontSize', 16, 'FontWeight', 'bold');
    
    %% Calculate aHR and aROHR
    dp_dCA = diff(avgPressure) / resolution; 
    dV_dCA = diff(volume') / resolution; 
    % Calculate aROHR
    aROHR_result = get_aROHR(avgPressure, volume', gamma);
    % Calculate aHR
    aHR_result = get_aHR(aROHR_result); 
 
    % Store the results
    aHR_all{rowIdx} = aHR_result; % Store aHR for each file
    aROHR_all{rowIdx} = aROHR_result; %Store aROHR for each file
    % Plot results
    %plot(Ca(:, 1), aROHR_result , 'LineWidth', 1.5, 'Color', colors(i, :), 'DisplayName', sprintf('(CA = %d°)', crankAngle));


    %% Save the Figure
    % Uncomment the following lines to save each figure as a PNG file
    % saveFileName = sprintf('Group_%s_Comparison_PV.png', uniqueID);
    % saveas(gcf, saveFileName);
    
end
colors = lines(height(T)); 
crank_angle_trimmed = crank_angle(1:end-1);
%% Loop to plot aHR
figure;
hold on;
for i = 1:length(aHR_all)
    uniqueID = T.UniqueID{i};
    plot(crank_angle_trimmed(:, 1), aHR_all{i}, 'LineWidth', 1.5, 'Color', colors(i, :), ...
        'DisplayName', sprintf(' %d', uniqueID));
end
xlabel('Crank Angle (°)');
ylabel('Apparent Heat Release [J]');
title('aHR for All Files');
legend('show');
grid on;
xlim([-45, 135]);
hold off;

%% Loop to plot aROHR
figure;
hold on;
crank_angles = 15:21;
for i = 1:length(aROHR_all)
    uniqueID = T.UniqueID{i};
    plot(crank_angle_trimmed(:, 1), aROHR_all{i}, 'LineWidth', 1.5, 'Color', colors(i, :), ...
        'DisplayName', sprintf(' %d', uniqueID));
end
xlabel('Crank Angle (°)');
ylabel('Apparent Heat Release [J]');
title('aROHR for All Files');
legend('show');
grid on;
xlim([-45, 135]);
hold off;

fprintf('Completed plotting for all groups.\n');


% ============================================
% End of Script
% ============================================
