function AFR = compute_AFR(y_CO2, y_CO, y_O2, y_N2)
    % Input:
    % y_CO2 - Volumetric fraction of CO2 (e.g., 0.10 for 10%)
    % y_CO  - Volumetric fraction of CO (e.g., 0.01 for 1%)
    % y_O2  - Volumetric fraction of O2 (e.g., 0.05 for 5%)
    % y_N2  - Volumetric fraction of N2 (should sum to 1 with other gases)
    %
    % Output:
    % fuel_mass_flow - Mass flow rate of fuel (kg/h)
    % air_mass_flow  - Mass flow rate of air (kg/h)
    % AFR            - Air-Fuel Ratio
    
    % Molar masses (g/mol)
    M_CO2 = 44.01;
    M_CO  = 28.01;
    M_O2  = 32.00;
    M_N2  = 28.01;
    
    % Molar mass of diesel fuel C12H23 (g/mol)
    M_fuel = (12 * 12.01) + (23 * 1.008);  % = 167.3 g/mol
    
    % Mass fractions of exhaust gases
    numerator = y_CO2 * M_CO2 + y_CO * M_CO + y_O2 * M_O2 + y_N2 * M_N2;
    w_CO2 = (y_CO2 * M_CO2) / numerator;
    w_CO  = (y_CO  * M_CO ) / numerator;
    w_O2  = (y_O2  * M_O2 ) / numerator;
    w_N2  = (y_N2  * M_N2 ) / numerator;
    
    % Assume total exhaust mass flow rate (kg/h)
    m_exhaust = 1;  % Arbitrary basis
    
    % Mass flow rates of exhaust components (kg/h)
    m_CO2 = w_CO2 * m_exhaust;
    m_CO  = w_CO  * m_exhaust;
    m_O2  = w_O2  * m_exhaust;
    m_N2  = w_N2  * m_exhaust;
    
    % Mass fractions of carbon in CO2 and CO
    w_C_in_CO2 = (12.01) / M_CO2;
    w_C_in_CO  = (12.01) / M_CO;
    
    % Total mass of carbon in exhaust (kg/h)
    m_C_exhaust = (m_CO2 * w_C_in_CO2) + (m_CO * w_C_in_CO);
    
    % Mass fraction of carbon in diesel fuel
    w_C_in_fuel = (12 * 12.01) / M_fuel;
    
    % Mass flow rate of fuel (kg/h)
    fuel_mass_flow = m_C_exhaust / w_C_in_fuel;
    
    % Mass fraction of hydrogen in diesel fuel
    w_H_in_fuel = (23 * 1.008) / M_fuel;
    
    % Mass flow rate of hydrogen in fuel (kg/h)
    m_H_fuel = fuel_mass_flow * w_H_in_fuel;
    
    % Mass fraction of hydrogen in water
    w_H_in_H2O = (2 * 1.008) / 18.015;
    
    % Mass flow rate of water vapor (kg/h)
    m_H2O = m_H_fuel / w_H_in_H2O;
    
    % Mass fractions of oxygen in exhaust gases
    w_O_in_CO2 = (2 * 16.00) / M_CO2;
    w_O_in_CO  = (16.00) / M_CO;
    w_O_in_O2  = 1;  % Oxygen is pure O2
    w_O_in_H2O = (16.00) / 18.015;
    
    % Mass of oxygen in exhaust gases (kg/h)
    m_O_exhaust = (m_CO2 * w_O_in_CO2) + (m_CO * w_O_in_CO) + (m_O2 * w_O_in_O2) + (m_H2O * w_O_in_H2O);
    
    % Mass of oxygen consumed (kg/h) - assuming negligible oxygen in fuel
    m_O_consumed = m_O_exhaust;
    
    % Mass fraction of oxygen in air
    w_O_in_air = 0.23;  % Approximately 23% by mass
    
    % Mass flow rate of air (kg/h)
    air_mass_flow = m_O_consumed / w_O_in_air;
    
    % Air-Fuel Ratio (AFR)
    AFR = air_mass_flow / fuel_mass_flow;
    
    % Display results
    fprintf('Air-Fuel Ratio (AFR): %.2f\n', AFR);
end
