
%% %% Everything Related to aHR and aROHR
%% Rate of changes, Pressure and Volume
% Crank angle change per data point
dCA = resolution;

% Pressure change per data point
dp = diff(p_filtered_avg);

% Pressure change per Crank angle
dp_dCA = dp/dCA;

% Volume change per data point
dV = diff(V_avg);
% Volume change per Crank angle
dV_dCA = dV/dCA;

%% Calculate and Plot aROHR for All HVO50 Files
disp(size(p));
% Initialization
figure;
hold on;
colors = lines(size(session1_Raw_files, 1)); 
% Loop through each file
for i = 1:size(session1_Raw_files, 1)
    % Extracting file information
    filePath = session1_Raw_files{i, 1};
    crankAngle = session1_Raw_files{i, 3};

    % Load and preprocess data
    data = table2array(readtable(filePath)); 
    Ca = reshape(data(:, 1), [], Ncycles);      
    %pressure
    p = reshape(data(:, 2), [], Ncycles) * bara;  
    p_filtered = zeros(size(p));
   
    for j = 1:Ncycles
        p_filtered(:, j) = SGFilter(p(:, j), polynomialOrder, frameLength, 0);
    end

    p_filtered_avg = mean(p_filtered, 2);
    dp_dCA = diff(p_filtered_avg) / resolution; 
    dV_dCA = diff(V_avg) / resolution; 
    % Calculate aROHR
    aROHR_result = aROHR(p_filtered_avg, V_avg, resolution, gamma, dp_dCA, dV_dCA);
    % Calculate aHR
    aHR_result = aHR(aROHR_result, resolution); 

    % Store the results
    aHR_all{i} = aHR_result; % Store aHR for each file

    % Plot results
    plot(Ca(:, 1), aROHR_result , 'LineWidth', 1.5, 'Color', colors(i, :), 'DisplayName', sprintf('(CA = %d°)', crankAngle));
end

% Finalize plot
xlabel('Crank Angle (°)');
ylabel('Apparent Rate of Heat Release [J/°]');
title('aROHR for All HVO50 Files');
legend('show');
grid on;
xlim([-45, 135]); 
hold off;

% Plot aHR for all files
figure;
hold on;
crank_angles = 15:21;
for i = 1:length(aHR_all)

    plot(Ca(:, 1), aHR_all{i}, 'LineWidth', 1.5, 'Color', colors(i, :), ...
        'DisplayName', sprintf('CA %d', crank_angles(i)));
end
xlabel('Crank Angle (°)');
ylabel('Apparent Heat Release [J]');
title('aHR for All Files');
legend('show');
grid on;
xlim([-45, 135]);
hold off;


%% Calculate Apparent Heat Release
aHR_avg = aHR(aROHR_avg, resolution);  % Assuming aHR function is already defined

% Plot the Apparent Heat Release
figure;
plot(Ca(:, 1), aHR_avg, 'LineWidth', 1.5);  % Ca(:,1) is the crank angle array
xlabel('Crank Angle (°)');
ylabel('Apparent Heat Release [J]');
title('Apparent Heat Release (Average)');
xlim([-45,135]);
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



% % Ideal Diesel Cycle
% Extract initial conditions from actual data
% Find index corresponding to Intake Valve Closure (IVC)
% IVC_angle = ValveEvents.CaIVC;  % Crank angle for IVC
% P1 = p_filtered_avg(idx_IVC);  % Pressure at IVC
% V1 = V_avg(idx_IVC);           % Volume at IVC
% Assume inlet temperature or estimate based on conditions
% <<<<<<< HEAD:ExampleDataSet/Simple.m
% T1 = 295 ;  % K assumed to be ambient
% =======
% T1 = 298.15;  % K (Atmospheric Conditions)

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