function mfr_air = CalculateMassFlowAir(O2_percent, true_mfr_fuel, AFR_stoich)
% Inputs:
% O2_percent : Vector of O2% values from exhaust gas data.
% true_mfr_fuel   : Reshaped matrix of fuel mass flow rates (g/s).
% AFR_stoich : Stoichiometric air-fuel ratio.
%
% Output:
% mfr_air    : Calculated mass flow rate of air (kg/s).

    % Calculate % Excess Air
    excess_air = (O2_percent / (20.9 - O2_percent)) * 100;

    % Calculate equivilent ratio (Air Excess Ratio)
    equivilent_ratio = 1 + excess_air / 100;
    
    % Calculate Mass flow rate of air
    mfr_air = equivilent_ratio * AFR_stoich * true_mfr_fuel;

end