function [dataIn, ExhaustData, Ca, p_filt, S_current, mfr_fuel, CO_percent_load, HC_ppm_load, NOx_ppm_load, CO2_percent_load, O2_percent_load, lambda_load] = loadingfromT(T, ID, bara)

rowIndex = find(strcmp(T.UniqueID, ID ));   % Find the row index for the desired ID
averageCycleData = T.AverageCycleData{rowIndex};    % Extract the AverageCycleData for that row
experimentDataCell = T.ExperimentData{rowIndex};    % Get the 1x1 cell of ExperimentData
experimentData = experimentDataCell{1};         % Extract the 36000x4 matrix inside
filteredPressureCell = T.FilteredPressure{rowIndex};    %Get the 1x1 cell of the Filteredpressure
filteredPressure = filteredPressureCell{1};


Ca = experimentData(1:3600, 1);    % 1st column contains the crank angles
p = averageCycleData.AvgPressure; % 2nd column contains the pressure 
mfr_fuel = averageCycleData.AvgMassFlow; % 3rd column contains the mass flow rate of fuel
S_current = averageCycleData.AvgCurrent;  % 4th column contians the sensor current
p_filt = filteredPressure(1:3600, 1); % 5th column contains the filtered pressure data
dataIn = [Ca, p, mfr_fuel, S_current, p_filt];

% % I'm not sure if lines 18-23 are neccesary tbh (Kata 1/1/2024)
% Ncycles = 100;
% % Reshape data into cycles
% Ca = reshape(dataIn(:, 1), [], Ncycles);      % Crank angle in degrees
% p_filt = reshape(dataIn(:, 5), [], Ncycles) * bara;  % Pressure in Pa (FILTERED!)
% S_current = reshape(dataIn(:, 4), [], Ncycles);  % Sensor current 
% mfr_fuel = reshape(dataIn(:, 3), [], Ncycles);  % Fuel mass flow

%% Load in the exhaust data:
exhaustDatainT = T.AdditionalData{rowIndex};    % Extract the Exhaust data for the blend
CO_percent_load = exhaustDatainT.CO;
HC_ppm_load = exhaustDatainT.HC;
NOx_ppm_load = exhaustDatainT.NOx;
CO2_percent_load = exhaustDatainT.CO2;
O2_percent_load = exhaustDatainT.O2;
lambda_load = exhaustDatainT.Lambda;

ExhaustData = [CO_percent_load, HC_ppm_load, NOx_ppm_load, CO2_percent_load, O2_percent_load, lambda_load];

end