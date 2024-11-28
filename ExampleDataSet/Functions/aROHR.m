function aROHR = aROHR(p, V, gamma, dV_dtheta, dp_dtheta)
    % Function to compute the Apparent Rate of Heat Release (aROHR)
    %
    % Inputs:
    %   p          - Pressure (Pa)
    %   V          - Cylinder Volume (m^3)
    %   gamma      - Ratio of specific heats
    %   dV_dtheta  - Derivative of volume w.r.t. crank angle (m^3/deg)
    %   dp_dtheta  - Derivative of pressure w.r.t. crank angle (Pa/deg)
    %
    % Output:
    %   aROHR      - Apparent Rate of Heat Release (J/deg)


    
    % Calculate the aROHR
    aROHR = (gamma / (gamma - 1)) * p .* dV_dtheta + (1 / (gamma - 1)) * V .* dp_dtheta;
end


