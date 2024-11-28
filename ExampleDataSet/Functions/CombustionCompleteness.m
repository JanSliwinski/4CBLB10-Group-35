function [combustion_efficiency] = CombustionCompleteness(CO_wt_percent, CO2_wt_percent)
    % Molar masses (g/mol)
    M_CO = 28.01;
    M_CO2 = 44.01;
    
    % Assume total mass is 100 g
    total_mass = 100;
    
    % Masses (g)
    mass_CO = CO_wt_percent * total_mass / 100;
    mass_CO2 = CO2_wt_percent * total_mass / 100;
    
    % Moles
    n_CO = mass_CO / M_CO;
    n_CO2 = mass_CO2 / M_CO2;
    
    % Combustion efficiency (%)
    combustion_efficiency = n_CO2 / (n_CO2 + n_CO) * 100;
end
