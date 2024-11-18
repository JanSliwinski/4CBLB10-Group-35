function [V] = CylinderVolume(Ca, Cyl)
    % Inputs:
    %   Ca  : Crank angle [degrees]
    %   Cyl : Struct with engine geometry
    % Outputs:
    %   V   : Cylinder volume [m^3]
    
    % Extract parameters from the struct
    B = Cyl.Bore;               % Bore [m]
    S = Cyl.Stroke;             % Stroke [m]
    r = S / 2;                  % Crank radius [m]
    l = Cyl.ConRod;             % Connecting rod length [m]
    cr = Cyl.CompressionRatio;  % Compression ratio
    A = pi / 4 * B^2;           % Cross-sectional area [m^2]
    
    % Volumes
    Vd = A * S;                 % Displacement volume [m^3]
    Vc = Vd / (cr - 1);         % Clearance volume [m^3]
    
    % Adjust crank angle to be relative to TDC
    CAl = Ca;
    
    % Calculate cylinder volume as a function of crank angle
    V = Vc + A * (r * (1 - cosd(CAl)) + l - sqrt(l^2 - (r * sind(CAl)).^2));
end
