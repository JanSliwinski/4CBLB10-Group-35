function mfr_air = CalculateMassFlowAir(O2_percent, mfr_fuel, AFR_stoich)
% CalculateMassFlowAir computes the air mass flow rate based on measured O2%
% in the exhaust, the fuel mass flow rate, and the stoichiometric AFR.
%
% Inputs:
%   O2_percent : Vector or matrix of measured O2% in exhaust gas data.
%   mfr_fuel   : Vector or matrix of fuel mass flow rates (g/s).
%   AFR_stoich : Stoichiometric air-fuel ratio (scalar or vector).
%
% Output:
%   mfr_air    : Calculated mass flow rate of air (g/s). Has the same
%                dimensions as the inputs if they are compatible.

    % Calculate % Excess Air (using element-wise division)
    excess_air = (O2_percent ./ (20.9 - O2_percent)) * 100;

    % Calculate the air excess ratio (sometimes called 'lambda')
    % Note: lambda = (actual AFR) / (stoichiometric AFR)
    % Here, 'equivalent_ratio' = 1 + (excess_air / 100).
    equivalent_ratio = 1 + (excess_air ./ 100);

    % Calculate Mass flow rate of air (element-wise multiplication)
    % mfr_air = lambda * AFR_stoich * mfr_fuel
    mfr_air = equivalent_ratio .* AFR_stoich .* mfr_fuel;
end