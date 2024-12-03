function [stoich_coeffs, reaction_eq, AFR_stoich] = StoichiometricCombustion(fuel_name, SpS, El)
    % Inputs:
%   - fuel_name: (string) The name of the fuel species (e.g., 'Methane').
%   - SpS: (struct array) A database containing species properties, including:
%       * SpS.Name: (string) Name of the species (e.g., 'O2', 'CO2').
%       * SpS.Elcomp: (array) Elemental composition of the species [C, H, O, ...].
%   - El: (struct array) A database containing element properties, including:
%       * El.Name: (string) Name of the element (e.g., 'C', 'H', 'O').

% Outputs:
%   - stoich_coeffs: (struct) A structure containing stoichiometric coefficients:
%       * stoich_coeffs.fuel: Coefficient for the fuel.
%       * stoich_coeffs.O2: Coefficient for oxygen in reactants.
%       * stoich_coeffs.N2: Coefficient for nitrogen in reactants.
%       * stoich_coeffs.CO2: Coefficient for carbon dioxide in products.
%       * stoich_coeffs.H2O: Coefficient for water in products.
%       * stoich_coeffs.N2_prod: Coefficient for nitrogen in products.
%   - reaction_eq: (string) The balanced stoichiometric combustion reaction as a string.

%% Explanation
% 1. **Identify the Fuel and Its Elemental Composition:**
%    - Locate the fuel in the SpS database by matching its name.
%    - Extract the elemental composition of the fuel (number of C, H, and O atoms).
%    - Raise an error if the fuel does not contain C or H.

% 2. **Determine Stoichiometric Combustion Coefficients:**
%    - For a fuel molecule CxHyOz, compute:
%      - b = x (moles of CO2 in products).
%      - c = y / 2 (moles of H2O in products).
%      - a = (2b + c - z) / 2 (moles of O2 required for combustion).
%    - Assume air consists of 3.76 moles of N2 for every mole of O2.
%      - N2 in reactants = 3.76 * a.
%      - N2 in products = N2 in reactants.

% 3. **Convert to Integer Coefficients:**
%    - Scale all coefficients to integers by finding the least common multiple (LCM).

% 4. **Generate the Balanced Reaction Equation:**
%    - Construct a string representation of the reaction equation using the 
%      computed stoichiometric coefficients.
    
% Find the fuel in SpS
    fuel_idx = find(strcmp({SpS.Name}, fuel_name));
    if isempty(fuel_idx)
        error('Fuel %s not found in SpS.', fuel_name);
    end
    fuel_species = SpS(fuel_idx);

    % Extract elemental composition
    elcomp = fuel_species.Elcomp;  % Array of elemental composition

    % Extract element names and create a mapping
    element_names = {El.Name};

    % Find indices of elements
    idx_C = find(strcmp(element_names, 'C'));
    idx_H = find(strcmp(element_names, 'H'));
    idx_O = find(strcmp(element_names, 'O'));

    % Handle cases where elements might not be present
    if isempty(idx_C) || isempty(idx_H)
        error('Fuel must contain both Carbon and Hydrogen.');
    end
    if isempty(idx_O)
        idx_O = 0;
    end

    % Extract number of atoms from elcomp
    x = elcomp(idx_C);  % Number of Carbon atoms
    y = elcomp(idx_H);  % Number of Hydrogen atoms
    if idx_O > 0
        z = elcomp(idx_O);  % Number of Oxygen atoms
    else
        z = 0;
    end

    % Continue with the stoichiometric calculations
    % Stoichiometric combustion equation:
    % CxHyOz + a O2 + a*3.76 N2 -> b CO2 + c H2O + a*3.76 N2

    b = x;          % CO2 coefficient
    c = y / 2;      % H2O coefficient

    % Total oxygen required on product side
    O_required = 2*b + c - z;
    a = O_required / 2;  % O2 coefficient

    % Multiply to eliminate fractions and get integer coefficients
    coeffs = [1, a, b, c];
    denom_lcm = lcm(denominator(coeffs(2)), lcm(denominator(coeffs(3)), denominator(coeffs(4))));
    if denom_lcm == 0
        multiplier = 1;
    else
        multiplier = denom_lcm;
    end

    % Apply multiplier to get integer coefficients
    fuel_coeff = multiplier;
    O2_coeff = multiplier * a;
    N2_coeff = O2_coeff * 3.76;
    CO2_coeff = multiplier * b;
    H2O_coeff = multiplier * c;
    N2_prod_coeff = N2_coeff;

    % Store coefficients
    stoich_coeffs.fuel = fuel_coeff;
    stoich_coeffs.O2 = O2_coeff;
    stoich_coeffs.N2 = N2_coeff;
    stoich_coeffs.CO2 = CO2_coeff;
    stoich_coeffs.H2O = H2O_coeff;
    stoich_coeffs.N2_prod = N2_prod_coeff;
    % Reaction equation string
    reaction_eq = sprintf('%g %s + %g O2 + %g N2 -> %g CO2 + %g H2O + %g N2', ...
        fuel_coeff, fuel_name, O2_coeff, N2_coeff, CO2_coeff, H2O_coeff, N2_prod_coeff);

    % Molecular weights
    MW_O2 = 32;  % g/mol
    MW_N2 = 28;  % g/mol
    MW_air = MW_O2 + 3.76 * MW_N2;  % Air molecular weight
    MW_fuel = x * 12 + y * 1 + z * 16;  % Fuel molecular weight

    mass_air = a * MW_air;  % Mass of air required

    % Stoichiometric air-fuel ratio
    AFR_stoich = mass_air / MW_fuel;
end

function denom = denominator(x)
    % Helper function to find the denominator of a fraction
    [~, d] = rat(x, 1e-12);
    denom = d;
end
