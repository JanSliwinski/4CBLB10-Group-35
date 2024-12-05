function mfr_air = CalculateMassFlowAir(O2_percent_vector, mfr_fuel, AFR_stoich)
% Inputs:
% O2_percent : Vector of O2% values from exhaust gas data.
% mfr_fuel   : Reshaped matrix of fuel mass flow rates (g/s).
% AFR_stoich : Stoichiometric air-fuel ratio.
%
% Output:
% mfr_air    : Calculated mass flow rate of air (kg/s).

    % Calculate % Excess Air
    excess_air = (O2_percent_vector ./ (20.9 - O2_percent_vector)) * 100;

    % Calculate equivilent ratio (Air Excess Ratio)
    equivilent_ratio = 1 + excess_air / 100;

    % Ensure Data Consistency
    if size(mfr_fuel, 2) ~= length(equivilent_ratio)
        error('Mismatch between the number of cycles in fuel data and O2 percentage data!');
    end

    % Initialize air mass flow rate matrix
    mfr_air = zeros(size(mfr_fuel));

    % Loop through each cycle
    for i = 1:size(mfr_fuel, 2)
        mfr_air(:, i) = equivilent_ratio(i) .* AFR_stoich .* mfr_fuel(:, i);
    end
end
