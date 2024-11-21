%% Initialization
clear; clc; close all;
addpath("Functions", "Nasa");  % Add necessary paths

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

%% Load NASA Data (if needed)
global Runiv
Runiv = 8.314;

[SpS, El] = myload('Nasa/NasaThermalDatabase.mat', {'Diesel', 'O2', 'N2', 'CO2', 'H2O'});

%% Engine Geometry Data
Cyl.Bore = 104 * mm;
Cyl.Stroke = 85 * mm;
Cyl.CompressionRatio = 21.5;
Cyl.ConRod = 136.5 * mm;
Cyl.TDCangle = 180;

%% Valve Events (Crank Angles)
ValveEvents.CaIVO = -355;
ValveEvents.CaIVC = -135;
ValveEvents.CaEVO = 149;
ValveEvents.CaEVC = -344;
ValveEvents.CaSOI = -3.2;  % Start of Injection

%% Load and Reshape Data
dataFileName = fullfile('Data', 'ExampleDataSet.txt');
dataIn = table2array(readtable(dataFileName));

[Nrows, ~] = size(dataIn);
resolution = 0.2;  % Degrees crank angle resolution
NdatapointsPerCycle = 720 / resolution;
Ncycles = Nrows / NdatapointsPerCycle;

if mod(Nrows, NdatapointsPerCycle) ~= 0
    error('Number of data points is not an integer multiple of data points per cycle.');
end

% Reshape data into cycles
Ca = reshape(dataIn(:, 1), [], Ncycles);      % Crank angle in degrees
p = reshape(dataIn(:, 2), [], Ncycles) * bara;  % Pressure in Pa

%% Plot Pressure vs. Crank Angle for All Cycles
figure;
set(gcf, 'Position', [200, 800, 1200, 400]);

plot(Ca, p / bara, 'LineWidth', 1);
xlabel('Crank Angle (°)');
ylabel('Pressure (bar)');
xlim([-360, 360]);
ylim([0, 50]);
title('Pressure vs. Crank Angle for All Cycles');
grid on;

% Highlight a specific cycle
iselect = 10;
hold on;
plot(Ca(:, iselect), p(:, iselect) / bara, 'r', 'LineWidth', 2);

% Plot valve events
YLIM = ylim;
line([ValveEvents.CaIVC, ValveEvents.CaIVC], YLIM, 'Color', 'b', 'LineWidth', 1);
line([ValveEvents.CaEVO, ValveEvents.CaEVO], YLIM, 'Color', 'r', 'LineWidth', 1);

set(gca, 'XTick', -360:60:360);
grid on;

%% Calculate Cylinder Volume for All Cycles
V_all = zeros(size(Ca));  % Initialize volume matrix
for i = 1:Ncycles
    V_all(:, i) = CylinderVolume(Ca(:, i), Cyl);
end

%% Filter Pressure Data
polynomialOrder = 3;
frameLength = 21;  % Must be odd

% Initialize the filtered pressure matrix
p_filtered = zeros(size(p));

% Apply the filter to each column
for i = 1:Ncycles
    p_filtered(:, i) = SGFilter(p(:, i), polynomialOrder, frameLength, 0);
end

%% Calculate Average Volume and Pressure
V_avg = mean(V_all, 2);         % Average volume across all cycles
p_avg = mean(p, 2);             % Average pressure across all cycles
p_filtered_avg = mean(p_filtered, 2);

%% Calculate Work
W = trapz(V_avg, p_avg);
disp(['Calculated work: ', num2str(W), ' J']);

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

%% Compare Raw and Filtered Pressure Data for a Single Cycle
figure;
hold on;
plot(Ca(:, iselect), p(:, iselect) / bara, 'DisplayName', 'Raw Data');
plot(Ca(:, iselect), p_filtered(:, iselect) / bara, 'DisplayName', 'Filtered Data');
xlabel('Crank Angle (°)');
ylabel('Pressure [bar]');
title(['Comparison of Raw and Filtered Pressure Data (Cycle ', num2str(iselect), ')']);
legend('show');
grid on;
hold off;

%% Parameters for the Ideal Diesel Cycle
% Extract initial conditions from actual data
% Find index corresponding to Intake Valve Closure (IVC)
IVC_angle = ValveEvents.CaIVC;  % Crank angle for IVC
[~, idx_IVC] = min(abs(Ca(:, 1) - IVC_angle));

P1 = p_filtered_avg(idx_IVC);  % Pressure at IVC
V1 = V_avg(idx_IVC);           % Volume at IVC

% Assume inlet temperature or estimate based on conditions
T1 = 300;  % K (adjust if you have actual data)

% Compression ratio
r = Cyl.CompressionRatio;

% Identify index of peak pressure in actual data
[~, idx_peak] = max(p_filtered_avg);

% Initial guesses and bounds for parameters
rc_initial_guess = 2.0;  % Increased initial guess for rc
T4_initial_guess = 1800;  % K
k_initial_guess = 1.35;

% Bounds for rc, T4, and k
lb = [1.1, 1000, 1.3];  % Lower bounds
ub = [3.0, 2500, 1.4];  % Upper bounds for rc increased to 3.0

numPoints = 100;  % Number of points per process

%% Optimization Function
totalErrorFunction = @(params) totalPressureDifference(params, Cyl, P1, T1, numPoints, V_avg, p_filtered_avg, idx_peak, idx_IVC);

% Optimize parameters
options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 1000);
optimized_params = fmincon(totalErrorFunction, [rc_initial_guess, T4_initial_guess, k_initial_guess], ...
                           [], [], [], [], lb, ub, [], options);

% Extract optimized parameters
rc = optimized_params(1);
T4 = optimized_params(2);
k = optimized_params(3);

% Display optimized parameters
disp('Optimized Parameters:');
disp(['Cutoff Ratio (rc): ', num2str(rc)]);
disp(['Maximum Temperature (T4): ', num2str(T4), ' K']);
disp(['Specific Heat Ratio (k): ', num2str(k)]);

% Recalculate the ideal cycle with optimized parameters
[P_cycle, V_cycle] = IdealDieselCycle(Cyl, P1, T1, T4, numPoints, k, rc);

% Convert volumes and pressures to match units in actual data
V_cycle_dm3 = V_cycle / dm^3;
P_cycle_bar = P_cycle / bara;

%% Plot the Adjusted Ideal Diesel Cycle Against Actual Data
% Plot filtered average pV diagram
figure;
loglog(V_avg / dm^3, p_filtered_avg / bara, 'b', 'LineWidth', 1.5);
hold on;

% Plot adjusted ideal Diesel cycle
loglog(V_cycle_dm3, P_cycle_bar, 'r--', 'LineWidth', 2);

% Plot key points of the ideal Diesel cycle
plot(V_cycle_dm3(1), P_cycle_bar(1), 'ko', 'MarkerFaceColor', 'k');                       % Point 1
plot(V_cycle_dm3(numPoints), P_cycle_bar(numPoints), 'go', 'MarkerFaceColor', 'g');       % Point 2
plot(V_cycle_dm3(2*numPoints), P_cycle_bar(2*numPoints), 'ro', 'MarkerFaceColor', 'r');   % Point 3
plot(V_cycle_dm3(3*numPoints), P_cycle_bar(3*numPoints), 'mo', 'MarkerFaceColor', 'm');   % Point 4

xlabel('Volume [dm³]');
ylabel('Pressure [bar]');
title('Optimized Ideal Diesel Cycle with Key Points');
legend('Filtered Average Data', 'Optimized Ideal Diesel Cycle', 'Point 1', 'Point 2', 'Point 3', 'Point 4');
grid on;
hold off;

%% Function Definitions

% Function to calculate the total pressure difference over the cycle
function error = totalPressureDifference(params, Cyl, P1, T1, numPoints, V_avg, p_filtered_avg, idx_peak, idx_IVC)
    rc = params(1);
    T4 = params(2);
    k = params(3);
    
    [P_cycle, V_cycle] = IdealDieselCycle(Cyl, P1, T1, T4, numPoints, k, rc);
    
    % Ensure V_avg and p_filtered_avg have unique, sorted volume values
    [V_avg_unique, idx_unique] = unique(V_avg, 'sorted');
    p_filtered_avg_unique = p_filtered_avg(idx_unique);
    
    % Interpolate actual pressure at theoretical volumes
    p_actual_interp = interp1(V_avg_unique, p_filtered_avg_unique, V_cycle, 'linear', 'extrap');
    
    % Calculate the error as the sum of squared differences
    pressure_error = sum((P_cycle - p_actual_interp).^2);
    
    % Calculate the volume change during combustion in the ideal cycle
    delta_V_combustion = (rc - 1) * V_cycle(numPoints);  % V3 - V2
    
    % Calculate the expected volume change from actual data
    delta_V_actual = V_avg(idx_peak) - V_avg(idx_IVC);
    
    % Calculate the difference in volume change
    volume_error = (delta_V_combustion - delta_V_actual)^2;
    
    % Combine errors with a weighting factor
    error = pressure_error + 1000 * volume_error;  % Adjust weighting factor as needed
end

% IdealDieselCycle Function (should be saved in a separate file IdealDieselCycle.m or included here)
function [P_cycle, V_cycle] = IdealDieselCycle(Cyl, P1, T1, T4, numPoints, k, rc)
    % Cylinder dimensions
    Bore = Cyl.Bore;
    Stroke = Cyl.Stroke;
    r = Cyl.CompressionRatio;  % Compression ratio
    
    % Calculate displacement volume Vd
    Vd = (pi / 4) * Bore^2 * Stroke;
    
    % Clearance volume Vc
    Vc = Vd / (r - 1);
    
    % Volumes at key points
    V1 = Vc + Vd;  % Volume at Bottom Dead Center (BDC)
    V2 = Vc;       % Volume at Top Dead Center (TDC)
    
    % Process 1-2: Isentropic compression (from V1 to V2)
    T2 = T1 * (V1 / V2)^(k - 1);
    P2 = P1 * (V1 / V2)^k;
    
    % Volume at point 3 (after heat addition)
    V3 = rc * V2;
    
    % Process 2-3: Constant-pressure heat addition (from V2 to V3)
    P3 = P2;  % Pressure remains constant
    T3 = T2 * (V3 / V2);  % T3 = T2 * (V3 / V2)
    
    % Process 3-4: Isentropic expansion (from V3 to V4, where V4 = V1)
    V4 = V1;
    T4_check = T3 * (V4 / V3)^(k - 1);  % Should match T4
    P4 = P3 * (V3 / V4)^k;
    
    % Verify T4
    if abs(T4 - T4_check) > 1e-3
        warning('Calculated T4 does not match the given T4.');
    end
    
    % Generate arrays for each process
    % Process 1-2: Isentropic compression
    V_12 = linspace(V1, V2, numPoints);
    P_12 = P1 * (V1 ./ V_12).^k;
    
    % Process 2-3: Constant-pressure heat addition
    V_23 = linspace(V2, V3, numPoints);
    P_23 = P2 * ones(size(V_23));
    
    % Process 3-4: Isentropic expansion
    V_34 = linspace(V3, V4, numPoints);
    P_34 = P3 * (V3 ./ V_34).^k;
    
    % Process 4-1: Constant-volume heat rejection (from P4 to P1)
    V_41 = V4 * ones(1, numPoints);
    P_41 = linspace(P4, P1, numPoints);
    
    % Concatenate the arrays to complete the cycle
    V_cycle = [V_12, V_23, V_34, V_41];
    P_cycle = [P_12, P_23, P_34, P_41];
end

% CylinderVolume Function (should be provided in your codebase)
% Ensure this function calculates the instantaneous cylinder volume based on crank angle and cylinder geometry

% SGFilter Function (should be provided in your codebase)
% Ensure this function applies the Savitzky-Golay filter to smooth the pressure data
