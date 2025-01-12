function [dataIn, ExhaustData, Ca, p_avg, S_current, mfr_fuel, CO_percent_load, HC_ppm_load, NOx_ppm_load, CO2_percent_load, O2_percent_load, lambda_load,  x_blend, LHV_blend, perc_blend, fuel_used, M_blend] = loadingfromT(T, ID, bara, x_GTL, x_diesel, x_HVO, LHV_GTL, LHV_HVO, LHV_diesel, M_GTL, M_HVO, M_diesel)

rowIndex = find(strcmp(T.UniqueID, ID ));   % Find the row index for the desired ID
averageCycleData = T.AverageCycleData{rowIndex};    % Extract the AverageCycleData for that row
experimentDataCell = T.ExperimentData{rowIndex};    % Get the 1x1 cell of ExperimentData
experimentData = experimentDataCell{1};         % Extract the 36000x4 matrix inside

Ca = experimentData(1:3600, 1);    % 1st column contains the crank angles
p_avg = averageCycleData.AvgPressure; % 2nd column contains the pressure 
mfr_fuel = averageCycleData.AvgMassFlow; % 3rd column contains the mass flow rate of fuel
S_current = averageCycleData.AvgCurrent;  % 4th column contians the sensor current
dataIn = [Ca, p_avg, mfr_fuel, S_current];

%% Load in the exhaust data:
exhaustDatainT = T.AdditionalData{rowIndex};    % Extract the Exhaust data for the blend
CO_percent_load = exhaustDatainT.CO;
HC_ppm_load = exhaustDatainT.HC;
NOx_ppm_load = exhaustDatainT.NOx;
CO2_percent_load = exhaustDatainT.CO2;
O2_percent_load = exhaustDatainT.O2;
lambda_load = exhaustDatainT.Lambda;

ExhaustData = [CO_percent_load, HC_ppm_load, NOx_ppm_load, CO2_percent_load, O2_percent_load, lambda_load];

%% Extract information from the ID automatically
% ID follows the format: L-Load%, I-IgnitionCrankAngle, C-Blend%, Fuel-FuelName

% Extract Fuel name and Blend percentage from the ID
fuel_pattern = '(?<=Fuel)[A-Za-z0-9]+';
blend_pattern = '(?<=C)\d+';

% Extract Fuel
fuel_match = regexp(ID, fuel_pattern, 'match', 'once');

% Extract Blend Percentage
blend_match = regexp(ID, blend_pattern, 'match', 'once');
blend_value = str2double(blend_match);

% Assign fuel_used based on C and fuel_match
if blend_value == 0
    fuel_used = 'diesel';
elseif contains(fuel_match, 'GTL', 'IgnoreCase', true)
    fuel_used = sprintf('GTL%d', blend_value);
elseif contains(fuel_match, 'HVO', 'IgnoreCase', true)
    fuel_used = sprintf('HVO%d', blend_value);
else
    error('Fuel type not recognized');
end

% Calculate Blend Percentage
perc_blend = blend_value / 100;

% Determine the blend type (x_blend) and LHV_blend based on fuel_used
if contains(fuel_used, 'GTL', 'IgnoreCase', true)
    x_blend = x_GTL;
    LHV_blend = LHV_GTL;
    M_blend = M_GTL;
elseif contains(fuel_used, 'diesel', 'IgnoreCase', true)
    x_blend = x_diesel;
    LHV_blend = LHV_diesel;
    M_blend = M_diesel;
elseif contains(fuel_used, 'HVO', 'IgnoreCase', true)
    x_blend = x_HVO;
    LHV_blend = LHV_HVO;
    M_blend = M_HVO;
else
    error('Fuel type not recognized');
end

% Output the variables
fprintf('Fuel used: %s\n', fuel_used);
fprintf('Blend percentage: %.2f\n', perc_blend);

end