%% Example calculations showcasing the use of introduced functions 
% for load 3.5 CA = 14 (default)

% Clear command window, workspace, and close all figures
clc, clear all, close all; warning off

%% Add necessary paths
% Set relative path to NASA database folder
relativepath_to_generalfolder = 'Nasa'; % Adjust if necessary
addpath(relativepath_to_generalfolder);

%% Load Nasa database
% Construct full path to thermal database and load it
TdataBase = fullfile('Nasa', 'NasaThermalDatabase');
load(TdataBase);

%% Find species
% Locate indexes of specific species in the database
iSp = myfind({Sp.Name}, {'O2', 'CO2', 'N2','H2O'}); % Find indexes of these species
SpS = Sp(iSp); % Create subset of species based on found indexes
NSp = length(SpS); % Number of species

%% Given variables
% Define mole fractions for each species
% fractions represent typical air composition
O2_frac = 0.1447;      % Oxygen fraction
CO2_frac = 0.0467;    % Carbon dioxide fracitons
N2_frac = 0.7808;      % Nitrogen fraction
% Calculate water vapor fraciton by subtraction
H2O_frac = 0.2095 - O2_frac - CO2_frac;
CO_frac = 0;      % Carbon monoxide fraction

% Combustion and thermal parameters
% These should also be read from a table/structure.
LHV = 43 * 1e6;           % Lower Heating Value in J/kg
m_dot_fuel = 0.0013;      % Mass flow rate of fuel (kg/s)
Q_dot = LHV * m_dot_fuel; % Heat transfer rate (W)
T_initial = 295.15;       % Initial temperature (K)
tolerance = 1e2;          % Acceptable error in heat transfer (W)
deltaT = 100;             % Initial guess for temperature change (K)
error = 1000;              % Initialize error to infinite

% Calculate Air-Fuel Ratio (AFR)
AFR = compute_AFR(CO2_frac,CO_frac,O2_frac,N2_frac);
m_dot_air = AFR * m_dot_fuel;
m_dot_tot = m_dot_air + m_dot_fuel; % Total mass flow rate

%% Convert mole fractions to mass fractions

% Collect mole fractions into an array
moleFractions = [O2_frac, CO2_frac, N2_frac, H2O_frac];

% Get molar masses of species
Mi = [SpS.Mass]; % Molar masses in kg/mol

% Calculate mass composition based on mole fractions and molar masses
massComposition = (moleFractions .* Mi) / sum(moleFractions .* Mi);

% Initialize counter for iteration tracking
counter = 0;
cpCum = []; % Array to potentially store cumulative cp values
TCum = [];  % Array to potentially store cumulative temperatures

%% Iterative calculation to find temperature change
% Iterate until calculated heat transfer is within tolerance of given heat transfer
while error > tolerance
    % Calculate new temperature
    T = T_initial + deltaT;
    
    % Compute average temperature for cp calculation
    T_avg = (T + T_initial) / 2;
    
    % Calculate specific heat capacity at constant pressure
    cp = compute_cp(T_avg, SpS, massComposition);
    
    % Calculate heat transfer based on current cp and temperature change
    Q_dot_calculated = m_dot_tot * cp * deltaT;
    
    % Compute the error between calculated and given heat transfer
    error = abs(Q_dot_calculated - Q_dot);
    
    % Increment temperature change and iteration counter
    deltaT = deltaT + 0.1;
    counter = counter + 1;

    % for plotting
    TCum = [TCum T_avg];
    cpCum = [cpCum cp];
end

% Display success message and results
disp("Great Success, cp = ");
disp(cp);

plot(TCum,cpCum)
xlabel("Temperature")
ylabel("C_p")
title("c_p for post-combustion mixture vs temperature")
subtitle("linear behaviour?")

%% Repeated display of results
disp("Great Success, cp = ");
disp(cp);
disp("Counter = ")
disp(counter)
disp("DeltaT = ")
disp(deltaT)

%% Find temperature where specific heat capacity becomes negative
% % Start at a high temperature
% T_test = 5000;
% cp_test = 5;
% 
% % Decrease temperature until cp becomes negative
% while cp_test > 0
%     cp_test = compute_cp(T_test,SpS,massComposition);
%     T_test = T_test + 1;
% end
% 
% % Compute cp at a temperature slightly below the negative cp point
% cp_show = compute_cp(T_test-2000,SpS,massComposition);
% 
% % Display results of cp investigation
% fprintf("c_p for T > %i is negative! = %f \n", T_test, cp_test);
% fprintf("c_p for %i is %f \n", T_test-2000, cp_show);

%%
% Load preprocessed data from a specific file
data = load("./Data/Processed_session1/filtered_averaged_data_3.5_IMEP.txt");

% Crop data from 180*5 to 180*15 indices (focus on specific crank angle range)
data = data(180*5:180*15,:);

% Ambient and gas constant parameters
T_amb = 295;   % Ambient temperature [K]
p_amb = 1e5;   % Ambient pressure [Pa]
R_air = 287;   % Gas constant for air
R_exhaust = 292;  % Gas constant for exhaust
rho = 1.202;   % Air density
mm = 0.001;    % Millimeter conversion

% Engine Geometry Parameters
Cyl.Bore = 104 * mm;               % Cylinder bore
Cyl.Stroke = 85 * mm;              % Cylinder stroke
Cyl.CompressionRatio = 21.5;       % Compression ratio
Cyl.ConRod = 136.5 * mm;           % Connecting rod length
Cyl.TDCangle = 180;                % Top Dead Center angle

% Calculate cylinder volume using CylinderVolume function
data(:,5) = CylinderVolume(data(:,1),Cyl);

% Find maximum pressure and combustion point
[maxP,combust] = max(data(:,2));

% Calculate initial air mass
m_air = data(1,2)*1e5 * data(1,5) / (R_air * T_amb);

% Temperature and gamma calculation loop
for dummy = 1:size(data,1)
    % Pre-combustion calculations
    if dummy <= combust
        % Calculate temperature using ideal gas law
        data(dummy,6) = data(dummy,2)*1e5 * data(dummy,5) / (R_air * m_air);
        
        % Calculate specific heat ratio (gamma)
        data(dummy,7) = compute_cp(data(dummy,6),SpS,massComposition) ...
            / (compute_cp(data(dummy,6),SpS,massComposition) - R_air);
    
    % Post-combustion calculations
    else
        % Calculate temperature using exhaust gas properties
        data(dummy,6) = data(dummy,2)*1e5 * data(dummy,5) / (R_exhaust * m_air);
        
        % Calculate specific heat ratio (gamma)
        data(dummy,7) = compute_cp(data(dummy,6),SpS,massComposition) ...
            / (compute_cp(data(dummy,6),SpS,massComposition) - R_exhaust);
    end
end

% First figure: Temperature and Gamma vs Crank Angle
figure;
yyaxis left
plot(data(:,1), data(:,6), '-b', 'LineWidth', 1);
ylabel('Temperature [K]');
xlabel('CA [deg]');
grid on;

yyaxis right
plot(data(:,1), data(:,7), 'r', 'LineWidth', 1);
ylabel('Gamma [-]');
title('Temperature and Gamma vs Crank Angle');
legend({'Temperature', 'Gamma'}, 'Location', 'best');

% Second figure: Gamma vs Temperature (pre and post combustion)
figure;
plot(data(1:combust,6), data(1:combust,7))
hold on
plot(data(combust+1:end,6), data(combust+1:end,7))
xlabel('Temperature [K]');
ylabel('Gamma [-]');
title("Gamma vs Temperature")
legend(["Pre-combustion" "Post-combustion"])