function [P_cycle, V_cycle] = IdealDieselCycle(Cyl)
%IDEALDIESELCYCLE Computes the ideal Diesel engine cycle pressure and volume curves.
%
% Inputs:
%   Cyl      : A structure containing cylinder parameters:
%              - Bore             : Cylinder bore [m]
%              - Stroke           : Cylinder stroke [m]
%              - CompressionRatio : Compression ratio (dimensionless)
%
% Outputs:
%   P_cycle : Array of pressures [Pa] along the ideal Diesel cycle
%   V_cycle : Array of volumes [m^3] along the ideal Diesel cycle
%
% This function calculates the pressure-volume (P-V) diagram of an ideal Diesel
% cycle based on given cylinder geometry using assumed thermodynamic processes:
% isentropic compression, constant-pressure heat addition, isentropic expansion,
% and constant-volume heat rejection.

    %% Constants and Initial Conditions
    P1 = 101325;      % Initial pressure [Pa]
    T1 = 298.15;      % Initial temperature [K]
    rc = 2;           % Ratio for end of heat addition (assumed value)
    k = 1.312562;     % Specific heat ratio (assumed constant)
    numPoints = 100;  % Number of points for each process curve

    %% Cylinder Dimensions and Calculated Volumes
    Bore = Cyl.Bore;
    Stroke = Cyl.Stroke;
    r = Cyl.CompressionRatio;  % Compression ratio

    % Calculate displacement volume Vd
    Vd = (pi / 4) * Bore^2 * Stroke;
    
    % Calculate clearance volume Vc using compression ratio
    Vc = Vd / (r - 1);
    
    % Volumes at key positions
    V1 = Vc + Vd;  % Volume at Bottom Dead Center (BDC)
    V2 = Vc;       % Volume at Top Dead Center (TDC)

    %% Process 1-2: Isentropic Compression (BDC to TDC)
    % Compute temperature and pressure at end of compression
    T2 = T1 * (V1 / V2)^(k - 1);
    P2 = P1 * (V1 / V2)^k;
    
    %% Process 2-3: Constant-Pressure Heat Addition
    % Expand volume at constant pressure from V2 to V3
    V3 = rc * V2;  % Volume after heat addition
    P3 = P2;         % Pressure remains constant during 2-3
    T3 = T2 * (V3 / V2);  % Temperature after heat addition
    
    %% Process 3-4: Isentropic Expansion
    % Isentropic expansion from V3 back to V1
    V4 = V1;                
    P4 = P3 * (V3 / V4)^k;  % Pressure at end of expansion

    %% Generate Process Curves
    % Process 1-2: Isentropic compression from V1 to V2
    V_12 = linspace(V1, V2, numPoints);
    P_12 = P1 * (V1 ./ V_12).^k;
    
    % Process 2-3: Constant-pressure heat addition from V2 to V3
    V_23 = linspace(V2, V3, numPoints);
    P_23 = P2 * ones(size(V_23));
    
    % Process 3-4: Isentropic expansion from V3 to V4
    V_34 = linspace(V3, V4, numPoints);
    P_34 = P3 * (V3 ./ V_34).^k;
    
    % Process 4-1: Constant-volume heat rejection from pressure P4 to P1 at V4
    V_41 = V4 * ones(1, numPoints);
    P_41 = linspace(P4, P1, numPoints);
    
    %% Concatenate Cycle Data
    % Combine all process segments to form complete cycle arrays
    V_cycle = [V_12, V_23, V_34, V_41];
    P_cycle = [P_12, P_23, P_34, P_41];
end
