clc, clear, close all

%% Initial Parameters
bore = 0.104;                     % Diameter of piston [m]
stroke = 0.085;                   % Full travel of piston [m]
conrodLength = 0.1365;            % Length of rod connecting piston head and crankshaft [m]
cRatio = 21.5;                    % Compression ratio [-]
displacement = 722e-6;            % The volume of the cylinder [m^3]
P0 = 101325;                      % Ambient Pressure [Pa]
r = stroke / 2;

% Global constants
global Runiv Pref Tref
Runiv = 8.314472;                 % Universal gas constant
Pref = 1.01235e5;                 % Reference pressure, 1 atm!
Tref = 298.15;                    % Reference Temperature [K]

% Diesel properties (B7 diesel - commonly available in Europe)
rho = 836.1;                      % Density [kg/m^3]
cn = 52.2;                        % Cetane ratio (ignitability) [-]
eta = 2.7638e-6;                  % Viscosity [m^2/s]

%% Ideal Diesel Cycle
% Extract initial conditions from actual data
% Find index corresponding to Intake Valve Closure (IVC)
IVC_angle = ValveEvents.CaIVC;  % Crank angle for IVC
P1 = p_filtered_avg(idx_IVC);  % Pressure at IVC
V1 = V_avg(idx_IVC);           % Volume at IVC
% Assume inlet temperature or estimate based on conditions
T1 = 298.15;  % K (Atmospheric Conditions)


r = Cyl.CompressionRatio;% Compression ratio
numPoints = 100;  % Number of points per process
rc = optimized_params(1);
T4 = optimized_params(2);%Exhaust Temprature
k = optimized_params(3);%Specific Heat Ratio


% Calculate the ideal cycle
[P_cycle, V_cycle] = IdealDieselCycle(Cyl, P1, T1, T4, numPoints, k, rc);

% Convert volumes and pressures to match units in actual data
V_cycle_dm3 = V_cycle / dm^3;
P_cycle_bar = P_cycle / bara;


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
