function aHR = get_aHR(aROHR)
% aHR Apparent Heat Release
% Calculates the apparent heat release by integrating the aROHR over crank angle.
%
% Inputs:
%   aROHR     - Apparent rate of heat release [J/deg]
%   resolution - Crank angle step size [degrees]
%
% Output:
%   aHR - Apparent heat release [J]
    resolution = 0.2;

    % Perform cumulative trapezoidal integration
    aHR = cumtrapz(aROHR) .* resolution; % Integration over aROHR
end

