%% HeatOfCombustion Function Documentation
% This function calculates the heat of combustion (ΔH_combustion) for a specified fuel
% using its stoichiometric combustion coefficients and standard enthalpies of formation.
% The calculation assumes standard conditions (T = 298.15 K).

% Function Signature:
% [deltaH_combustion] = HeatOfCombustion(fuel_name, SpS, stoich_coeffs)

% Inputs:
%   - fuel_name: (string) The name of the fuel species (e.g., 'Methane').
%   - SpS: (struct array) A database containing species properties, including:
%       * SpS.Name: (string) Name of the species.
%       * SpS.Ts: (array) Temperature range(s) for NASA polynomials.
%       * SpS.Pol: (matrix) Polynomial coefficients for enthalpy and entropy.
%   - stoich_coeffs: (struct) A structure containing stoichiometric coefficients, including:
%       * stoich_coeffs.fuel: Coefficient for the fuel.
%       * stoich_coeffs.O2: Coefficient for oxygen.
%       * stoich_coeffs.N2: Coefficient for nitrogen in reactants.
%       * stoich_coeffs.CO2: Coefficient for carbon dioxide in products.
%       * stoich_coeffs.H2O: Coefficient for water in products.
%       * stoich_coeffs.N2_prod: Coefficient for nitrogen in products.

% Outputs:
%   - deltaH_combustion: (double) The heat of combustion in Joules (J).
%       * Negative value indicates an exothermic reaction.

%% Explanation
% 1. **Identify Species:**
%    - Locate the required species (O2, N2, CO2, H2O) and the fuel in the SpS database.
%    - Ensure all species are found, otherwise raise an error.

% 2. **Calculate Enthalpies at Standard Temperature (298.15 K):**
%    - Use the `SpeciesEnthalpy` helper function to calculate enthalpy for each species:
%      - Fuel
%      - O2 (reactant)
%      - N2 (reactant and product)
%      - CO2 (product)
%      - H2O (product)

% 3. **Reactants and Products Enthalpy:**
%    - Compute the total enthalpy of reactants
%      - Compute the total enthalpy of products

% 4. **Heat of Combustion:**
%    - Calculate the heat of combustion

function [deltaH_combustion] = HeatOfCombustion(fuel_name, SpS, stoich_coeffs)
    % Get species indices
    species_names = {'O2', 'N2', 'CO2', 'H2O'};
    species_indices = zeros(length(species_names), 1);
    for i = 1:length(species_names)
        idx = find(strcmp({SpS.Name}, species_names{i}));
        if isempty(idx)
            error('Species %s not found in SpS.', species_names{i});
        end
        species_indices(i) = idx;
    end

    % Get fuel index
    fuel_idx = find(strcmp({SpS.Name}, fuel_name));
    if isempty(fuel_idx)
        error('Fuel %s not found in SpS.', fuel_name);
    end

    T_std = 298.15;  % Standard temperature

    % Enthalpies of formation at standard temperature
    H_fuel = SpeciesEnthalpy(T_std, SpS(fuel_idx));
    H_O2 = SpeciesEnthalpy(T_std, SpS(species_indices(1)));
    H_N2 = SpeciesEnthalpy(T_std, SpS(species_indices(2)));
    H_CO2 = SpeciesEnthalpy(T_std, SpS(species_indices(3)));
    H_H2O = SpeciesEnthalpy(T_std, SpS(species_indices(4)));

    % Reactants and products enthalpy
    H_reactants = stoich_coeffs.fuel * H_fuel + stoich_coeffs.O2 * H_O2 + stoich_coeffs.N2 * H_N2;
    H_products = stoich_coeffs.CO2 * H_CO2 + stoich_coeffs.H2O * H_H2O + stoich_coeffs.N2_prod * H_N2;

    deltaH_combustion = H_products - H_reactants;  % Should be negative (exothermic)
end

function h = SpeciesEnthalpy(T, species)
% This helper function calculates the enthalpy of a species at a given temperature
% using NASA polynomials.

% Function Signature:
% h = SpeciesEnthalpy(T, species)

% Inputs:
%   - T: (double) Temperature in Kelvin (K).
%   - species: (struct) A structure representing the species, containing:
%       * species.Name: (string) Name of the species.
%       * species.Ts: (array) Temperature range(s) for NASA polynomials.
%       * species.Pol: (matrix) Polynomial coefficients for enthalpy.
%       * species.Hf: (double, optional) Standard enthalpy of formation (J/mol).

% Outputs:
%   - h: (double) Enthalpy in Joules (J/mol).
%% Explanation
% 1. **Temperature Range Selection:**
%    - Determine whether the input temperature T falls in the low or high range.
%    - Select the corresponding set of NASA polynomial coefficients.

% 2. **Enthalpy Calculation:**
%    - Use the NASA polynomial equation:
%    - Scale temperature t = T / 1000 to match NASA polynomial conventions.

% 3. **Error Handling:**
%    - Raise an error if T  is out of the valid range for the species.
    
%    - Use the standard enthalpy of formation H_f for species at 298.15 K when available.
    R = 8.314;  % J/mol·K

    T_ranges = species.Ts;  % Temperature ranges
    Pol = species.Pol;      % Polynomial coefficients

    % Check if species has multiple temperature ranges
    if length(T_ranges) == 3
        T_low = T_ranges(1);
        T_common = T_ranges(2);
        T_high = T_ranges(3);
    else
        % Assume only one temperature range
        T_low = T_ranges(1);
        T_common = T_ranges(1);
        T_high = T_ranges(end);
    end

    % Determine which set of coefficients to use
    if T >= T_low && T <= T_high
        if T < T_common && size(Pol, 1) >= 2
            a = Pol(2, :);  % Low-temperature coefficients
        else
            a = Pol(1, :);  % High-temperature coefficients
        end

        t = T / 1000;  % Scale temperature

        % Enthalpy calculation using NASA polynomials
        h_RT = a(1) + a(2)*t/2 + a(3)*t^2/3 + a(4)*t^3/4 + a(5)*t^4/5 + a(6)/t;
        h = h_RT * R * T;
    else
        % If temperature is outside the range, use standard enthalpy of formation if available
        if isfield(species, 'Hf') && T == 298.15
            h = species.Hf;
        elseif strcmp(species.Name, 'Diesel') && T == 298.15
            h = -249e3;  % Standard enthalpy of formation in J/mol
        else
            error('Temperature %.2f K is out of range for species %s.', T, species.Name);
        end
    end
end

