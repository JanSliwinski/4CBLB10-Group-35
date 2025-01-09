%function gamma = CalculateGamma(T,expNo)
clear all; close all; clc; expNo = 3; load("T.mat");
Runiv = 8.314; Tref = 293.15; pref = 1e5; mm=0.001;

%% Unpacking the table T
CA = (0:0.2:720-0.2)';
pressure = T.AverageCycleData{expNo}.AvgPressure;

%% Load Nasa database
addpath("Nasa\");
TdataBase = fullfile('Nasa', 'NasaThermalDatabase'); 
load(TdataBase);

%% Find species
iSp = myfind({Sp.Name}, {'O2', 'CO2', 'N2','H2O'}); % Find indexes of these species
SpS = Sp(iSp); % Create subset of species based on found indexes
NSp = length(SpS); % Number of species

%% Exhuast composition data
N2_volfrac = 0.78;          % Nitrogen fraction in air
O2_volfrac_air = 0.21;      % Oxygen fraction in air
O2_volfrac = T.AdditionalData{expNo}.O2 / 100; % Oxygen fraction
CO2_volfrac = T.AdditionalData{expNo}.CO2 / 100; % Carbon dioxide fracitons
CO_volfrac = T.AdditionalData{expNo}.CO / 100; % Carbon monoxide fraction
H2O_volfrac = O2_volfrac_air - O2_volfrac - CO2_volfrac; % Calculate water vapor fraciton by subtraction

%% Compute specific gas constants for inlet and outlet
% Get molar masses of species in kg/mol
Mi = [SpS.Mass];
% Find specific gas constants for all constituents
R_co2 = Runiv/Mi(2); R_o2 = Runiv/Mi(1); R_n2 = Runiv/Mi(3); R_h20 = Runiv/Mi(4);

% Find mass fractions
totalMass = N2_volfrac*Mi(3)+O2_volfrac*Mi(1)+H2O_volfrac*Mi(4)+CO2_volfrac*Mi(2);
N2_massfrac = N2_volfrac*Mi(3) / totalMass;
O2_massfrac = O2_volfrac*Mi(1) / totalMass;
CO2_massfrac = CO2_volfrac*Mi(2) / totalMass;
H2O_massfrac = H2O_volfrac*Mi(4) / totalMass;
massFractions = [O2_massfrac,CO2_massfrac,N2_massfrac,H2O_massfrac];

% Get specific gas constats
R_in = 287; % Ambient air
R_out = N2_massfrac*R_n2+O2_massfrac*R_o2+CO2_massfrac*R_co2+H2O_massfrac*R_h20;

%% getting air mass per cycle and cylinder volume
% Engine Geometry Parameters
Cyl.Bore = 104 * mm;               % Cylinder bore
Cyl.Stroke = 85 * mm;              % Cylinder stroke
Cyl.CompressionRatio = 21.5;       % Compression ratio
Cyl.ConRod = 136.5 * mm;           % Connecting rod length
Cyl.TDCangle = 180;                % Top Dead Center angle

% Calculate cylinder volume using CylinderVolume function
volume = CylinderVolume(CA,Cyl);

[maxVolume,idx] = max(volume);
massAir = pressure(idx)*1e5 * maxVolume / (R_in * Tref);

%% Getting Temperature

% Find maximum pressure and combustion point
[maxP,combust] = max(pressure);

temperature = zeros(size(pressure));
gamma = zeros(size(pressure));
% Temperature and gamma calculation loop
for dummy = 1:size(pressure)
    % Pre-combustion calculations
    if dummy <= combust
        % Calculate temperature using ideal gas law
        temperature(dummy) = pressure(dummy)*1e5 * volume(dummy) / (R_in * massAir);
        
        % Calculate specific heat ratio (gamma)
        gamma(dummy) = compute_cp(temperature(dummy),SpS,massFractions) ...
            / (compute_cp(temperature(dummy),SpS,massFractions) - R_in);
    
    % Post-combustion calculations
    else
        % Calculate temperature using exhaust gas properties
        temperature(dummy) = pressure(dummy)*1e5 * volume(dummy) / (R_out * massAir);
        
        % Calculate specific heat ratio (gamma)
        gamma(dummy) = compute_cp(temperature(dummy),SpS,massFractions) ...
            / (compute_cp(temperature(dummy),SpS,massFractions) - R_out);
    end
end
%% cHECKIGN
plot(CA,pressure.*volume.^gamma);

%% Plotting
% First figure: Temperature and Gamma vs Crank Angle
figure;
yyaxis left
plot(CA, temperature, '-b', 'LineWidth', 1);
ylabel('Temperature [K]');
xlabel('CA [deg]');
grid on;

yyaxis right
plot(CA, volume, 'r', 'LineWidth', 1);
ylabel('Gamma [-]');
title('Temperature and Gamma vs Crank Angle');
legend({'Temperature', 'Gamma'}, 'Location', 'best');

% Second figure: Gamma vs Temperature (pre and post combustion)
figure;
plot(temperature(1:combust), gamma(1:combust))
hold on
plot(temperature(combust+1:end), gamma(combust+1:end))
xlabel('Temperature [K]');
ylabel('Gamma [-]');
title("Gamma vs Temperature")
legend(["Pre-combustion" "Post-combustion"])


%end