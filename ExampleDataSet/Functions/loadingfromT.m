function [dataIn, ExhaustData, Ca, p_avg, S_current, mfr_fuel, CO_percent_load, HC_ppm_load, NOx_ppm_load, CO2_percent_load, O2_percent_load, lambda_load, x_blend, LHV_blend, perc_blend, fuel_used, M_blend] = loadingfromT(T, ID, bara, x_GTL, x_diesel, x_HVO, LHV_GTL, LHV_HVO, LHV_diesel, M_GTL, M_HVO, M_diesel)
%LOADINGFROMT Loads experimental and exhaust data given an identifier.
%
% Inputs:
%   T           : Data table containing experimental results and metadata
%   ID          : Unique identifier string for the desired experiment
%   bara        : Pressure conversion factor (bar to Pa)
%   x_GTL,x_diesel,x_HVO: Carbon atoms in respective fuels
%   LHV_GTL,LHV_HVO,LHV_diesel: Lower heating values (J/g) for respective fuels
%   M_GTL,M_HVO,M_diesel: Molar masses (g/mol) for respective fuels
%
% Outputs:
%   dataIn          : Matrix combining crank angle, average pressure, mass flow, sensor current
%   ExhaustData     : Array containing exhaust measurements [CO%, HC ppm, NOx ppm, CO2%, O2%, Lambda]
%   Ca              : Crank angle array
%   p_avg           : Averaged pressure array
%   S_current       : Sensor current array
%   mfr_fuel        : Averaged fuel mass flow rate
%   CO_percent_load : CO percentage in exhaust
%   HC_ppm_load     : HC concentration in exhaust [ppm]
%   NOx_ppm_load    : NOx concentration in exhaust [ppm]
%   CO2_percent_load: CO2 percentage in exhaust
%   O2_percent_load : O2 percentage in exhaust
%   lambda_load     : Lambda value from exhaust measurement
%   x_blend         : Carbon atom count in fuel blend
%   LHV_blend       : Lower heating value of the fuel blend (J/g)
%   perc_blend      : Blend fraction (decimal)
%   fuel_used       : Identifier for fuel type used
%   M_blend         : Molar mass of the fuel blend (g/mol)
%
% This function extracts experimental and exhaust data corresponding to a given ID 
% from table T, parses the ID to determine the fuel type and blend percentage, 
% and computes additional fuel properties based on the fuel type and blend.

    %% Extract Data Row
    % Find row corresponding to given ID in table T
    rowIndex = find(strcmp(T.UniqueID, ID));   
    
    %% Retrieve Experimental Data
    % Extract average cycle and raw experiment data matrices for the specified row
    averageCycleData = T.AverageCycleData{rowIndex};
    experimentDataCell = T.ExperimentData{rowIndex};
    experimentData = experimentDataCell{1};  % Extract the 36000x4 matrix
    
    %% Extract Signals from Experimental Data
    Ca = experimentData(1:3600, 1);                 % Crank angles [deg]
    p_avg = averageCycleData.AvgPressure;           % Averaged pressure [Pa]
    mfr_fuel = averageCycleData.AvgMassFlow;        % Averaged fuel mass flow [kg/s]
    S_current = averageCycleData.AvgCurrent;        % Averaged sensor current [A]
    dataIn = [Ca, p_avg, mfr_fuel, S_current];      % Combine signals for input
    
    %% Load Exhaust Data
    % Extract exhaust measurements from the table
    exhaustDatainT = T.AdditionalData{rowIndex};
    CO_percent_load = exhaustDatainT.CO;
    HC_ppm_load = exhaustDatainT.HC;
    NOx_ppm_load = exhaustDatainT.NOx;
    CO2_percent_load = exhaustDatainT.CO2;
    O2_percent_load = exhaustDatainT.O2;
    lambda_load = exhaustDatainT.Lambda;
    
    ExhaustData = [CO_percent_load, HC_ppm_load, NOx_ppm_load, CO2_percent_load, O2_percent_load, lambda_load];
    
    %% Parse ID to Determine Fuel Type and Blend
    % Extract fuel type and blend percentage from ID using regular expressions
    fuel_pattern = '(?<=Fuel)[A-Za-z0-9]+';  % Pattern for fuel name
    blend_pattern = '(?<=C)\d+';            % Pattern for blend percentage
    
    fuel_match = regexp(ID, fuel_pattern, 'match', 'once');    % Extract fuel type
    blend_match = regexp(ID, blend_pattern, 'match', 'once');  % Extract blend value
    blend_value = str2double(blend_match);                    % Convert blend string to number
    
    %% Determine Fuel Used Based on Blend Value and Fuel Match
    if blend_value == 0
        fuel_used = 'diesel';
    elseif contains(fuel_match, 'GTL', 'IgnoreCase', true)
        fuel_used = sprintf('GTL%d', blend_value);
    elseif contains(fuel_match, 'HVO', 'IgnoreCase', true)
        fuel_used = sprintf('HVO%d', blend_value);
    else
        error('Fuel type not recognized');
    end
    
    %% Calculate Blend Percent and Fuel Properties
    perc_blend = blend_value / 100;  % Convert blend percentage to fraction
    
    % Determine fuel properties based on identified fuel type
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
    
    %% Output Fuel Information
    fprintf('Fuel used: %s\n', fuel_used);
    fprintf('Blend percentage: %.2f\n', perc_blend);
end
