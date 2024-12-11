function aROHR = aROHR(p_filtered_avg, V_avg, resolution, gamma, dp_dCA, dV_dCA)
% aROHR Apparent Rate of Heat Release
% Calculates the apparent rate of heat release using pressure and volume data.
%
% Inputs:
%   p_filtered_avg - Filtered average pressure [Pa]
%   V_avg          - Average volume [m³]
%   resolution     - Crank angle step size [degrees]
%   gamma          - Specific heat ratio (assumed constant)
%   dp_dCA         - Pressure direvative over Crank angle [Pa/degree]
%   dV_dCA         - Volume direvative over Crank angle [m³/degree]
%
% Output:
%   aROHR - Apparent rate of heat release [J/deg]

    % Convert resolution to radians
    % dTheta = deg2rad(resolution);

    % Differentiate pressure with respect to crank angle
    % dp_dTheta = diff(p_filtered_avg) ./ dTheta;

    % Differentiate volume with respect to crank angle
    % dV_dTheta = diff(V_avg) ./ dTheta;

    % Calculate ROHR using the formula
    aROHR = gamma/(gamma - 1) * p_filtered_avg(1:end-1) .* dV_dCA + (1 / (gamma - 1)) * (V_avg(1:end-1) .* dp_dCA);

    % Append zero for the last point to match dimensions
    aROHR = [aROHR; 0];
end
