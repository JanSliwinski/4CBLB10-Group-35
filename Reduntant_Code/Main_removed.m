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