function aHR = aHR(aROHR, resolution)
% aHR Apparent Heat Release
% Calculates the apparent heat release by integrating the aROHR over crank angle.
%
% Inputs:
%   aROHR     - Apparent rate of heat release [J/deg]
%   resolution - Crank angle step size [degrees]
%
% Output:
%   aHR - Apparent heat release [J]

    % Convert resolution to radians for consistency
    dTheta = deg2rad(resolution);

    % Perform cumulative trapezoidal integration
    aHR = cumtrapz(aROHR) * dTheta; % Integration over aROHR
end

