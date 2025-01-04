function aROHR = get_aROHR(p_filt, volume, gamma)
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
    
        pressure = p_filt .* 1e5;
        
        resolution = 0.2;
        dp_dCA = [diff(pressure) / resolution;0]; 
        dV_dCA = [diff(volume) / resolution;0]; 
        aROHR = gamma./(gamma - 1) .* pressure .* dV_dCA + (1 ./ (gamma - 1)) .* (volume .* dp_dCA);
end
