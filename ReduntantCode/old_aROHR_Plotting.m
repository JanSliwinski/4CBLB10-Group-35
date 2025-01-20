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