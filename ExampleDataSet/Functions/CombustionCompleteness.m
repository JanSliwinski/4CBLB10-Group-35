function [combustion_efficiency] = CombustionCompleteness(CO_vol_percent, CO2_vol_percent)
    % Combustion completeness based on volume percent
    % Input:
    %   CO_vol_percent: Volume percent of CO (0-100)
    %   CO2_vol_percent: Volume percent of CO2 (0-100)
    % Output:
    %   combustion_efficiency: Combustion efficiency as a percentage (0-100)

    % Convert volume percent to volume fraction
    vol_fraction_CO = CO_vol_percent / 100;
    vol_fraction_CO2 = CO2_vol_percent / 100;

    % Calculate combustion efficiency
    combustion_efficiency = vol_fraction_CO2 / (vol_fraction_CO2 + vol_fraction_CO) * 100;
end
