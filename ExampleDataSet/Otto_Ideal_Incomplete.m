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

%% Assumed Data
lhv = [43e6; 37e6];               % Lower heating value [J/kg] (left: pure diesel, right: FAME)
RPMn = 1000;                      % Unloaded RPM of a typical diesel engine
rotTime = 60 / RPMn;              % Time of a full rotation (360 degrees) [s]
R_specific = 287;                 % Specific gas constant for air [J/(kg·K)]
md = 0.00005;                     % Fuel mass per cycle [kg]
Cp = 2.75e3;                      % Specific heat capacity at high temperature [J/kg·K]
gamma = 1.4;                      % Gamma for air-fuel mixture
Cv = Cp / gamma;                  % Specific heat capacity at constant volume
ma = 14.5 * md;                   % Air-fuel mass ratio (14.5:1)
mtot = ma + md;                   % Total mass (air + fuel)

% Crank angle and time setup
t = linspace(0, rotTime, 720);    % Time vector [s]
ca = 1:720;                       % Crank angle in degrees

% Initialize pressure, temperature, and volume
p = zeros(1, 720); T = zeros(1, 720); V = zeros(1, 720);
p(1) = P0; T(1) = Tref;           % Initial pressure and temperature
V(1) = Volume(ca(1), cRatio, conrodLength, bore, displacement, r);

% Heat release during combustion
dQ_comb = sum(lhv .* [md * 0.93; md * 0.07]); % Heat from diesel (93%) and FAME (7%)

% Loop over each crank angle
for n = 2:720
    V(n) = Volume(ca(n), cRatio, conrodLength, bore, displacement, r);

    switch true
        %% Intake Stroke
        case (ca(n) <= 180)
            p(n) = P0;
            T(n) = Tref;

        %% Compression Stroke
        case (ca(n) > 180 && ca(n) <= 360)
            T(n) = T(n-1) * (V(n-1) / V(n))^(gamma - 1);
            p(n) = p(n-1) * (V(n-1) / V(n))^gamma;

        %% Combustion Phase (Isobaric for Diesel Cycle)
        case (ca(n) > 360 && ca(n) <= 380)
            p(n) = p(n-1); % Isobaric combustion
            T(n) = T(n-1) + (dQ_comb / (mtot * Cv));

        %% Expansion Stroke
        case (ca(n) > 380 && ca(n) <= 540)
            T(n) = T(n-1) * (V(n-1) / V(n))^(gamma - 1);
            p(n) = p(n-1) * (V(n-1) / V(n))^gamma;

        %% Exhaust Phase
        case (ca(n) > 540 && ca(n) <= 720)
            p(n) = P0;
            T(n) = Tref;

        otherwise
            disp('Unexpected crank angle detected!');
    end
end

%% Plotting Results

% P-V Diagram
figure(1)
plot(V(2:end), p(2:end), 'LineWidth', 2);
grid on;
xlabel('Volume (m^3)');
ylabel('Pressure (Pa)');
title('P-V Diagram');

% Temperature vs. Crank Angle
figure(2)
plot(ca, T, 'LineWidth', 2);
grid on;
xlabel('Crank Angle (degrees)');
ylabel('Temperature (K)');
title('Temperature vs. Crank Angle');

% Log-log P-V Diagram
figure(3)
loglog(V(2:end), p(2:end), 'LineWidth', 2);
grid on;
xlabel('log(Volume)');
ylabel('log(Pressure)');
title('Log-Log P-V Diagram');

%% Volume Function
function V = Volume(ca, cRatio, conrodLength, bore, displacement, r)
    % Input parameters:
    % ca: Crank angle (in degrees)
    % cRatio: Compression ratio
    % conrodLength: Length of the connecting rod (m)
    % bore: Bore diameter of the engine (m)
    % displacement: Total displacement volume (m^3)

    % Convert crank angle to radians
    ca = deg2rad(ca);

    % Clearance volume based on compression ratio
    Vc = displacement / (cRatio - 1);

    % Volume calculation at a given crank angle
    V = Vc + (pi / 4) * bore^2 * (r + conrodLength - (r * cos(ca) + sqrt(conrodLength^2 - r^2 * sin(ca)^2)));
end
