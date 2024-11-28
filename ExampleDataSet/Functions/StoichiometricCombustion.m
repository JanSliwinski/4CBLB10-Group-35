function [stoich_coeffs, reaction_eq] = StoichiometricCombustion(fuel_name, SpS, El)
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
end

function denom = denominator(x)
    % Helper function to find the denominator of a fraction
    [~, d] = rat(x, 1e-12);
    denom = d;
end
