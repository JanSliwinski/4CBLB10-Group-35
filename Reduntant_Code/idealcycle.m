clc;
clear;
close all;

% Constants and Inputs
gamma = 1.312562;                % Specific heat ratio (air)
R = 287;                    % Specific gas constant for air [J/kg·K]
LHV = 43e6;                 % Lower Heating Value of diesel [J/kg]
bore = 0.104;               % Bore diameter [m]
stroke = 0.085;             % Stroke length [m]
CR = 21.5;                  % Compression ratio
p1 = 101325;                % Initial pressure [Pa]
T1 = 298.15;                % Initial temperature [K]

% Geometric Calculations
V_swept = pi * (bore / 2)^2 * stroke;         
V_clearance = V_swept / (CR - 1);            
V1 = V_swept + V_clearance;                   % Volume at the start of compression [m^3]
V2 = V_clearance;                             % Volume at the end of compression [m^3]

% State 2: Adiabatic Compression
p2 = p1 * (V1 / V2)^gamma;                    % Pressure after compression [Pa]
T2 = T1 * (V1 / V2)^(gamma - 1);              % Temperature after compression [K]

% State 3: Isobaric Heat Addition
r_c = 2;                                     
V3 = r_c * V2;                                % Volume after combustion [m^3]
cp = R * gamma / (gamma - 1);                
T3 = T2 * r_c;                                
p3 = p2;                                      % Pressure remains constant during isobaric process

% State 4: Adiabatic Expansion
V4 = V1;                                      % Volume returns to maximum at end of expansion
p4 = p3 * (V3 / V4)^gamma;                    % Pressure after expansion [Pa]
T4 = T3 * (V3 / V4)^(gamma - 1);              % Temperature after expansion [K]

% Generate Points for Each Process
num_points = 100;

% Compression: From V1 to V2
V_comp = linspace(V1, V2, num_points);
P_comp = p1 * (V1 ./ V_comp).^gamma;

% Heat Addition (Isobaric): From V2 to V3
V_add = linspace(V2, V3, num_points);
P_add = repmat(p3, 1, num_points);

% Expansion: From V3 to V4
V_exp = linspace(V3, V4, num_points);
P_exp = p3 * (V3 ./ V_exp).^gamma;

% Heat Rejection (Isochoric): From V4 to V1
V_reject = repmat(V1, 1, num_points);
P_reject = linspace(p4, p1, num_points);

% Combine Data for Plotting
V_cycle = [V_comp, V_add, V_exp, V_reject];
P_cycle = [P_comp, P_add, P_exp, P_reject];

% Plot Normal PV Diagram (Linear scale)
figure;
plot(V_cycle, P_cycle / 1e6, 'LineWidth', 2);
hold on;
scatter([V1, V2, V3, V4], [p1, p2, p3, p4] / 1e6, 'o', 'filled');
text(V1, p1 / 1e6, '1', 'VerticalAlignment', 'bottom');
text(V2, p2 / 1e6, '2', 'VerticalAlignment', 'bottom');
text(V3, p3 / 1e6, '3', 'VerticalAlignment', 'bottom');
text(V4, p4 / 1e6, '4', 'VerticalAlignment', 'bottom');
xlabel('Volume [m³]');
ylabel('Pressure [MPa]');
title('Ideal Diesel Cycle PV Diagram (Linear Scale)');
grid on;
hold off;

% Plot Log-Log PV Diagram
figure;
loglog(V_cycle, P_cycle / 1e6, 'LineWidth', 2);  
hold on;
scatter([V1, V2, V3, V4], [p1, p2, p3, p4] / 1e6, 'o', 'filled');
text(V1, p1 / 1e6, '1', 'VerticalAlignment', 'bottom');
text(V2, p2 / 1e6, '2', 'VerticalAlignment', 'bottom');
text(V3, p3 / 1e6, '3', 'VerticalAlignment', 'bottom');
text(V4, p4 / 1e6, '4', 'VerticalAlignment', 'bottom');
xlabel('Volume [m³]');
ylabel('Pressure [MPa]');
title('Ideal Diesel Cycle PV Diagram (Log-Log Scale)');
grid on;
hold off;

% Diesel Cycle Efficiency
efficiency = 1 - (1 / CR^(gamma - 1)) * ((r_c^gamma - 1) / (gamma * (r_c - 1)));
fprintf('Thermal Efficiency: %.2f%%\n', efficiency * 100);
