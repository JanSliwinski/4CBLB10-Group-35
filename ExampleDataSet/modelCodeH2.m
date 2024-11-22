clear all;close all;clc;
warning off
%% To make sure that matlab will find the functions. You must change it to your situation 
relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase);
%% Nasa polynomials are loaded and globals are set. 
%% values should not be changed. These are used by all Nasa Functions. 
global Runiv Pref Tref
Runiv=8.314472;
Pref=1.01235e5; % Reference pressure, 1 atm!
Tref=298.15;    % Reference Temperature
%% Some convenient units
kJ=1e3;kmol=1e3;dm=0.1;bars=1e5;kPa = 1000;kN=1000;kg=1;s=1;
%% Given conditions. 
%  For the final assignment take the ones from the specific case you are supposed to do.                  
v1=200;Tamb=300;P3overP2=9.00;Pamb=100*kPa;mfurate=0.58*kg/s;AF=204.42;             % These are the ones from the book
cFuel='H2';                                                           % Pick Gasoline as the fuel (other choices check Sp.Name)
%% Select species for the case at hand
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
%% Air composition
Xair = [0 0.21 0 0 0.79];                                                   % Order is important. Note that these are molefractions
Mair = Xair*Mi';                                                            % Row times Column = inner product 
Yair = Xair.*Mi/Mair;                                                       % Vector. times vector is Matlab's way of making an elementwise multiplication
%% Fuel composition
Yfuel = [1 0 0 0 0];                                                        % Only fuel
%% Range of enthalpies/thermal part of entropy of species
TR = 200:1:3000;NTR=length(TR);
for i=1:NSp                                                                 % Compute properties for all species for temperature range TR 
    hia(:,i) = HNasa(TR,SpS(i));                                            % hia is a NTR by 5 matrix
    sia(:,i) = SNasa(TR,SpS(i));                                            % sia is a NTR by 5 matrix
end
hair_a= Yair*hia';                                                          % Matlab 'inner product': 1x5 times 5xNTR matrix muliplication, 1xNTR resulT -> enthalpy of air for range of T 
sair_a= Yair*sia';                                                          % same but this thermal part of entropy of air for range of T
% whos hia sia hair_a sair_a                                                  % Shows dimensions of arrays on commandline
%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the interpolation method
% Bisection is in the next 'cell'
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using INTERPOLATION
cMethod = 'Interpolation Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/Mair;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
T2 = interp1(hair_a,TR,h2);                                                 % Interpolate h2 on h2air_a to approximate T2. Pretty accurate
%knowing the h2, they are tryng to find the T2 basedon the avg enthalpy of
%air mixture
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
h2check = Yair*hi2';                                                        % Single value (1x5 times 5x1). Why do I do compute this h2check value? Any ideas?
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral part of th eentropy)
Pr = exp(lnPr);
P2 = P1*Pr;
s1  = s1thermal - Rg*log(P1/Pref);                                          % Total specific entropy
s2  = s2thermal - Rg*log(P2/Pref);
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',s1/kJ,s2/kJ);
T2int = T2;

%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the Bisection method
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using bisection (https://en.wikipedia.org/wiki/Bisection_method)
cMethod = 'Bisection Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/Mair;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
TL = T1;
TH = 1000;                                                                  % A guess for the TH (must be too high)
iter = 0;
while abs(TH-TL) > 0.01
    iter = iter+1;
    Ti = (TL+TH)/2;
    for i=1:NSp
        hi2(i)    = HNasa(Ti,SpS(i));
    end
    h2i = Yair*hi2';                                                        % Single value (1x5 times 5x1). Intermediate value
    if h2i > h2
        TH = Ti; % new right boundary
    else
        TL = Ti; % new left boundary
    end
end
T2 = (TH+TL)/2;
T2bis = T2;
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral)
Pr = exp(lnPr);
P2 = P1*Pr;
s1  = s1thermal - Rg*log(P1/Pref);                                          % Total entropy stage 1
s2  = s2thermal - Rg*log(P2/Pref);                                          % Total entropy stage 2
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',s1/kJ,s2/kJ);
%% Difference between two approaches: so close but not identical
fprintf('----------------------------------------------\n%8s| %9.4f %9.4f  [K]\n----------------------------------------------\n','T2-int vs T2-bis',T2int,T2bis);
%% Here starts your part (compressor,combustor,turbine and nozzle). ...
% Make a choice for which type of solution method you want to use.

%% Compressor 
Compressor = 'Compressor';
P3 = P3overP2 * P2;  % Calculate P3 using pressure ratio
v3 = v2;  % Velocity assumed to be zero in the jet engine untill nozzle
s3 = s2;  % Isentropic compression

s3thermal = s2thermal + Rg * log(P3overP2);  % Thermal entropy change

T3 = interp1(sair_a, TR, s3thermal);  % Interpolate temperature at stage 3

for i = 1:NSp
    hi3(i) = HNasa(T3, SpS(i));  % Calculate enthalpy at T3 for each species
end

h3 = Yair * hi3';  % Calculate total enthalpy at stage 3

% Print stage information for Compressor
fprintf('\n%14s\n','');
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',Compressor,2,3);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T2,T3);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P2/kPa,P3/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v2,v3);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h2/kJ,h3/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',s2/kJ,s3/kJ);

%% Combustor
Combustor = 'Combustor';
v4 = v3;
P4 = P3;  % As per the assignment's assumptions

mf_air = AF * mfurate;   % Mass flow rate of air
mf_total = mf_air + mfurate; % Total mass flow rate into the combustor

molmass_H2 = Mi(1);  % Molar mass of hydrogen (H2)
molmass_O2 = Mi(2);  % Molar mass of oxygen (O2)
molmass_N2 = Mi(5);  % Molar mass of nitrogen (N2)
molmass_H2O = Mi(4); % Molar mass of water (H2O)

% Molar flow rates entering the system
molflowin_H2 = mfurate / molmass_H2;                   % Molar flow rate in of hydrogen (H2)
molflowin_O2 = mf_air * Yair(1,2) / molmass_O2;        % Molar flow rate in of oxygen (O2)
molflowin_N2 = mf_air * Yair(1,5) / molmass_N2;        % Molar flow rate in of nitrogen (N2)

% Molar flow rates exiting the system
molflowout_H2O = molflowin_H2;                          % Molar flow rate out of water (H2O)
molflowout_O2 = molflowin_O2 - 0.5 *  molflowout_H2O;   % Molar flow rate out of oxygen (O2)
molflowout_N2 = molflowin_N2;                           % Molar flow rate out of nitrogen (N2, remains unchanged)

% Total molar flow rates
molflowin_tot = molflowin_H2 + molflowin_O2 + molflowin_N2; % Total molar flow rate in
molflowout_tot = molflowout_H2O + molflowout_O2 + molflowout_N2; % Total molar flow rate out

% Mass flow rates
mfin_H2 = mfurate;                                 % Mass flow rate in of hydrogen (H2)
mfin_O2 = mf_air * Yair(1,2);                      % Mass flow rate in of oxygen (O2)
mfin_N2 = mf_air * Yair(1,5);                      % Mass flow rate in of nitrogen (N2)

mfout_H2O = molflowout_H2O * molmass_H2O;           % Mass flow rate out of water (H2O)
mfout_O2 = molflowout_O2 * molmass_O2;              % Mass flow rate out of oxygen (O2)
mfout_N2 = molflowout_N2 * molmass_N2;              % Mass flow rate out of nitrogen (N2)


% Molar and mass fractions before and after combustion
Y3 = [mfin_H2 mfin_O2 0 0 mfin_N2]/mf_total;  % Mass fractions before combustion
X4 = [0 molflowout_O2/molflowout_tot 0 molflowout_H2O/molflowout_tot molflowout_N2/molflowout_tot]; % Mole fractions after combustion
Y4 = [0 mfout_O2 0 mfout_H2O mfout_N2]/mf_total;  % Mass fractions after combustion

% Getting elemental composition of fuel, oxygen, and water
fuel_comp = [SpS(1).Elcomp];    % Composition of fuel (H2)
O2_comp = [SpS(2).Elcomp];      % Composition of oxygen (O2)
H2O_comp = [SpS(4).Elcomp];     % Composition of water (H2O)

Mair_after = X4 * Mi';         % Molar mass of air after combustion
Rg_after = Runiv / Mair_after; % Gas constant after combustion

% Stoichiometric air-fuel ratio calculations
comp_ma_O2 = (fuel_comp(2) * 0.5) / 2;            % O2 moles needed for 1 mole of H2
comp_ma_N2 = Xair(5) / Xair(2) * comp_ma_O2;      % N2 moles needed for 1 mole of H2
stoiAFR = (comp_ma_N2 * Mi(5) + comp_ma_O2 * Mi(2)) / Mi(1);  % Stoichiometric air-fuel ratio (AFR)

% Equivalence ratio
equivalence_ratio = AF / stoiAFR;  % Equivalence ratio (actual AFR / stoichiometric AFR)

% Calculate the internal energy of the mixture before combustion
for i = 1:NSp
    ui1(i) = UNasa(T3, SpS(i));  % Internal energy of each species at T3
end

u_mix1 = Yair * ui1';  % Total internal energy of the mixture before combustion

% Compute internal energy after combustion for each temperature in TR
for i = 1:NSp
    ui2(:, i) = UNasa(TR, SpS(i));  % Internal energy for each species at various temperatures
end

u_mix2 = Y4 * ui2';  % Total internal energy of the mixture after combustion

% Interpolate to find T4 based on internal energy conservation
T4 = interp1(u_mix2, TR, u_mix1);

% Compute entropy and enthalpy for each species at T4
for i = 1:NSp
    si4(i) = SNasa(T4, SpS(i));  % Entropy of each species at T4
    hi4(i) = HNasa(T4, SpS(i));  % Enthalpy of each species at T4
end

h4 = Y4 * hi4';  % Total enthalpy at stage 4

s4thermal = Y4 * si4';  % Total thermal entropy

s4 = s4thermal - s3thermal + s3;  % Update total entropy

% Output stage information for Combustor
fprintf('\n%14s\n','');
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n', Combustor, 3, 4);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp', T3, T4);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press', P3/kPa, P4/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v', v3, v4);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h', h3/kJ, h4/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S', s3/kJ, s4/kJ);

%% Turbine
Turbine = 'Turbine';

v5 = v4;
s5 = s4;  % isentropic expansion

h5 = h4 + h2 - h3;  % Calculate enthalpy at stage 5 using energy balance

% Range of enthalpies across TR for the mixture after combustion 
h4_mix = Y4 * hia';

T5 = interp1(h4_mix, TR, h5);  % Interpolate temperature at stage 5

% Get the specific entropies for each species
for i = 1 : NSp
    si5(i) = SNasa(T5, SpS(i));
end

% Compute the total thermal entropy at stage 5
s5thermal = Y4 * si5';

% Find pressure based on thermal entropies
P5 = P4 * exp((s5thermal - s4thermal) / Rg_after);

% Output stage information for Turbine
fprintf('\n%14s\n','');
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n', Turbine, 4, 5);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n', 'Temp', T4, T5);
fprintf('%8s| %9.2f %9.2f  [kPa]\n', 'Press', P4/kPa, P5/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n', 'v', v4, v5);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n', 'h', h4/kJ, h5/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n', 'Total S', s4/kJ, s5/kJ);

%% Nozzle
Nozzle = 'Nozzle';

P6 = Pamb;  % Set exit pressure to ambient
s6 = s5;    % Entropy at stage 6 initialized as equal to stage 5

s6thermal = s5thermal + Rg_after * log(P6 / P5);  % Thermal entropy change in the nozzle

% Range of entropies across TR for the mixture after combustion
s4_mix = Y4 * sia'; 

T6 = interp1(s4_mix, TR, s6thermal);  % Interpolate temperature at stage 6

% Calculate enthalpy at stage 6 for each species
for i = 1:NSp
    hi6(i) = HNasa(T6, SpS(i));
end

h6 = Y4 * hi6';  % Total enthalpy at stage 6

v6 = sqrt(2 * (h5 - h6));  % Calculate exit velocity using enthalpy difference

% Output stage information for Nozzle
fprintf('\n%14s\n','');
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n', Nozzle, 5, 6);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n', 'Temp', T5, T6);
fprintf('%8s| %9.2f %9.2f  [kPa]\n', 'Press', P5/kPa, P6/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n', 'v', v5, v6);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n', 'h', h5/kJ, h6/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n', 'Total S', s5/kJ, s6/kJ);

%% Graphing
% Visualising the results, so we can easier grasp the jet engine
% Define state points for entropy vs temperature
TS1 = [s1 T1];
TS2 = [s2 T2];
TS3 = [s3 T3];
TS4 = [s4 T4];
TS5 = [s5 T5];
TS6 = [s6 T6];

entropies = [s1 s2 s3 s4 s5 s6];
temperatures = [T1 T2 T3 T4 T5 T6];

% Define state points for volume vs pressure
V1 = Rg * T1 / P1;
V2 = Rg * T2 / P2;
V3 = Rg * T3 / P3;
V4 = Rg * T4 / P4;
V5 = Rg * T5 / P5;
V6 = Rg * T6 / P6;

pressures = [P1 P2 P3 P4 P5 P6];
volumes = [V1 V2 V3 V4 V5 V6];

% Create a figure for subplots
figure

% First subplot: Entropy vs Temperature
subplot(1, 2, 1)
scatter(entropies, temperatures, 20, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')  % Red filled circles
text(entropies + 0.01, temperatures, {'1', '2', '3', '4', '5', '6'}, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')  % State labels
xlabel('Entropy [kJ/kgK]')
ylabel('Temperature [K]')
title('Entropy vs Temperature')
legend('States 1-6')
grid on

% Second subplot: Volume vs Pressure
subplot(1, 2, 2)
scatter(volumes, pressures, 20, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')  % Red filled circles
text(volumes + 0.01, pressures, {'1', '2', '3', '4', '5', '6'}, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')  % State labels
xlabel('Specific Volume [m^3/kg]')
ylabel('Pressure [kPa]')
title('Pressure vs Specific Volume')
legend('States 1-6')
grid on

% Overall title for the figure
sgtitle('Thermodynamic State Representations')

% Obtaining the power generatoin
P = 0.5* mf_total * v6^2;

% Obtaining thrust
F = P / v6;