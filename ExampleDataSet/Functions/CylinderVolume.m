function [V] = CylinderVolume(Ca, Cyl)
    % This function calculates the cylinder volume as a function of 
    % crank angle (Ca) and cylinder geometry.
    %
    % Inputs:
    %   Ca : Crankangle [degrees]
    %   Cyl : A struct containing:
    %       Cyl.S  : Stroke
    %       Cyl.B  : Bore
    %       Cyl.ConRod : Connecting Rod length
    %       Cyl.CompressionRatio : Compression Ratio
    %       Cyl.TDCangle : Angle associated with the Top Dead Center
    %
    % Outputs:
    %   V : Cylinder volume [cubic units]
    %----------------------------------------------------------------------
    
    B   = Cyl.Bore;
    S   = Cyl.S;
    cr  = Cyl.CompressionRatio;
    r   = S / 2; % Crank radius
    l   = Cyl.ConRod; % Connecting rod length
    
    % Calculate minimum and maximum cylinder volumes
    Vdiff = pi / 4 * B^2 * S; % Displacement volume
    Vmin = Vdiff / (cr - 1); % Clearance volume
    
    % Calculate instantaneous volume
    x = r + l - (r * cosd(Ca) + sqrt(l^2 - r^2 * (sind(Ca)).^2));
    V = Vmin + x * pi / 4 * B^2;
end
