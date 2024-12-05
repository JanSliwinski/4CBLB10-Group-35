%--------------------------------------------------------------------------
% Function to Calculate Mass Flow Rates of Fuel and Air Based on Exhaust Gas Composition
% and Determine the Air-Fuel Ratio (AFR) for a Diesel Engine Using Diesel Fuel (C12H23)
%--------------------------------------------------------------------------

function AFR = compute_AFR(y_CO2, y_CO, y_O2, y_N2)
    % Inputs:
    %   y_CO2 - Volumetric fraction of CO2 in the exhaust gas (e.g., 0.10 for 10%)
    %   y_CO  - Volumetric fraction of CO in the exhaust gas (e.g., 0.01 for 1%)
    %   y_O2  - Volumetric fraction of O2 in the exhaust gas (e.g., 0.05 for 5%)
    %   y_N2  - Volumetric fraction of N2 in the exhaust gas (should sum to 1 with other gases)
    %
    % Outputs:
    %   fuel_mass_flow - Calculated mass flow rate of fuel (kg/h)
    %   air_mass_flow  - Calculated mass flow rate of air (kg/h)
    %   AFR            - Air-Fuel Ratio (dimensionless)
    
    %--------------------------------------------------------------------------
    % Step 1: Define Molar Masses of Exhaust Gas Components and Fuel
    %--------------------------------------------------------------------------
    % Molar masses are in grams per mole (g/mol)
    
    % Molar masses of exhaust gas components
    M_CO2 = 44.01;    % Carbon dioxide (CO2)
    M_CO  = 28.01;    % Carbon monoxide (CO)
    M_O2  = 32.00;    % Oxygen (O2)
    M_N2  = 28.01;    % Nitrogen (N2)
    
    % Molar mass of diesel fuel approximated as C12H23
    % Calculated by summing the molar masses of carbon and hydrogen atoms
    M_fuel = (12 * 12.01) + (23 * 1.008);  % Diesel fuel (C12H23)
    % M_fuel = 167.3 g/mol
    
    %--------------------------------------------------------------------------
    % Step 2: Convert Volumetric Fractions to Mass Fractions
    %--------------------------------------------------------------------------
    % Under ideal gas conditions, volumetric fractions are equivalent to molar fractions.
    % To work with mass flow rates, we need to convert these to mass fractions using molar masses.
    
    % Calculate the weighted sum of molar masses multiplied by their respective volumetric fractions
    numerator = y_CO2 * M_CO2 + y_CO * M_CO + y_O2 * M_O2 + y_N2 * M_N2;
    
    % Compute mass fractions for each exhaust gas component
    w_CO2 = (y_CO2 * M_CO2) / numerator;  % Mass fraction of CO2
    w_CO  = (y_CO  * M_CO ) / numerator;  % Mass fraction of CO
    w_O2  = (y_O2  * M_O2 ) / numerator;  % Mass fraction of O2
    w_N2  = (y_N2  * M_N2 ) / numerator;  % Mass fraction of N2
    
    %--------------------------------------------------------------------------
    % Step 3: Assume Total Exhaust Mass Flow Rate (kg/h)
    %--------------------------------------------------------------------------
    % We assume an arbitrary total exhaust mass flow rate to establish a basis for calculation.
    % This simplifies the calculations and allows us to determine relative mass flow rates.
    
    m_exhaust = 1;  % Assumed total exhaust mass flow rate (kg/h)
    
    %--------------------------------------------------------------------------
    % Step 4: Calculate Mass Flow Rates of Exhaust Components
    %--------------------------------------------------------------------------
    % Multiply the mass fractions by the assumed total exhaust mass flow rate
    % to obtain mass flow rates for each component.
    
    m_CO2 = w_CO2 * m_exhaust;  % Mass flow rate of CO2 (kg/h)
    m_CO  = w_CO  * m_exhaust;  % Mass flow rate of CO (kg/h)
    m_O2  = w_O2  * m_exhaust;  % Mass flow rate of O2 (kg/h)
    m_N2  = w_N2  * m_exhaust;  % Mass flow rate of N2 (kg/h)
    
    %--------------------------------------------------------------------------
    % Step 5: Calculate Mass of Carbon in Exhaust Gases
    %--------------------------------------------------------------------------
    % Determine the mass fractions of carbon in CO2 and CO and calculate the
    % total mass of carbon leaving in the exhaust gases.
    
    % Mass fractions of carbon in CO2 and CO
    w_C_in_CO2 = 12.01 / M_CO2;  % Carbon mass fraction in CO2
    w_C_in_CO  = 12.01 / M_CO;   % Carbon mass fraction in CO
    
    % Total mass of carbon in the exhaust gases (kg/h)
    m_C_exhaust = (m_CO2 * w_C_in_CO2) + (m_CO * w_C_in_CO);
    
    %--------------------------------------------------------------------------
    % Step 6: Calculate Mass Flow Rate of Fuel Based on Carbon Mass Balance
    %--------------------------------------------------------------------------
    % Use the conservation of mass for carbon to determine the fuel mass flow rate.
    % All carbon in the exhaust is assumed to come from the fuel.
    
    % Mass fraction of carbon in diesel fuel (C12H23)
    w_C_in_fuel = (12 * 12.01) / M_fuel;
    
    % Fuel mass flow rate (kg/h)
    fuel_mass_flow = m_C_exhaust / w_C_in_fuel;
    
    %--------------------------------------------------------------------------
    % Step 7: Calculate Mass Flow Rate of Hydrogen in Fuel
    %--------------------------------------------------------------------------
    % Determine the mass flow rate of hydrogen in the fuel using the mass fraction.
    
    % Mass fraction of hydrogen in diesel fuel (C12H23)
    w_H_in_fuel = (23 * 1.008) / M_fuel;
    
    % Mass flow rate of hydrogen in fuel (kg/h)
    m_H_fuel = fuel_mass_flow * w_H_in_fuel;
    
    %--------------------------------------------------------------------------
    % Step 8: Calculate Mass Flow Rate of Water Vapor Produced
    %--------------------------------------------------------------------------
    % Assuming all hydrogen in the fuel is converted to water vapor in the exhaust.
    
    % Molar mass of water (H2O)
    M_H2O = 18.015;  % g/mol
    
    % Mass fraction of hydrogen in water
    w_H_in_H2O = (2 * 1.008) / M_H2O;
    
    % Mass flow rate of water vapor produced (kg/h)
    m_H2O = m_H_fuel / w_H_in_H2O;
    
    %--------------------------------------------------------------------------
    % Step 9: Calculate Mass of Oxygen in Exhaust Gases
    %--------------------------------------------------------------------------
    % Calculate the mass fractions of oxygen in each exhaust component and
    % compute the total mass of oxygen in the exhaust gases.
    
    % Mass fractions of oxygen in exhaust gases
    w_O_in_CO2 = (2 * 16.00) / M_CO2;  % Oxygen mass fraction in CO2
    w_O_in_CO  = 16.00 / M_CO;         % Oxygen mass fraction in CO
    w_O_in_O2  = 1.0;                  % Oxygen mass fraction in O2 (pure oxygen)
    w_O_in_H2O = 16.00 / M_H2O;        % Oxygen mass fraction in water
    
    % Total mass of oxygen in exhaust gases (kg/h)
    m_O_exhaust = (m_CO2 * w_O_in_CO2) + ...
                  (m_CO  * w_O_in_CO ) + ...
                  (m_O2  * w_O_in_O2 ) + ...
                  (m_H2O * w_O_in_H2O);
    
    %--------------------------------------------------------------------------
    % Step 10: Calculate Mass of Oxygen Consumed
    %--------------------------------------------------------------------------
    % Assuming negligible oxygen content in the fuel, all oxygen in the exhaust
    % gases is considered to come from the air consumed in combustion.
    
    m_O_consumed = m_O_exhaust;  % Mass of oxygen consumed (kg/h)
    
    %--------------------------------------------------------------------------
    % Step 11: Calculate Mass Flow Rate of Air Based on Oxygen Consumed
    %--------------------------------------------------------------------------
    % Use the mass fraction of oxygen in air to calculate the total air mass flow rate.
    
    % Mass fraction of oxygen in air (approximately 23% by mass)
    w_O_in_air = 0.23;
    
    % Air mass flow rate (kg/h)
    air_mass_flow = m_O_consumed / w_O_in_air;
    
    %--------------------------------------------------------------------------
    % Step 12: Calculate Air-Fuel Ratio (AFR)
    %--------------------------------------------------------------------------
    % The Air-Fuel Ratio is the ratio of the mass flow rate of air to the mass
    % flow rate of fuel.
    
    AFR = air_mass_flow / fuel_mass_flow;
    
    %--------------------------------------------------------------------------
    % Step 13: Display Results
    %--------------------------------------------------------------------------
    % Output the calculated mass flow rates and Air-Fuel Ratio.

    fprintf('Air-Fuel Ratio (AFR): %.2f\n', AFR);
end
