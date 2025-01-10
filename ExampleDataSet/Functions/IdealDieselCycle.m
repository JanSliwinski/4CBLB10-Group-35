  function [P_cycle, V_cycle] = IdealDieselCycle(Cyl)
    % Constants
    % k is provided as an input
    P1 = 101325;                % Initial pressure [Pa]
    T1 = 298.15;                % Initial temperature [K]
    rc = 2
    k = 1.312562
    numPoints = 100
    % Cylinder dimensions
    Bore = Cyl.Bore;
    Stroke = Cyl.Stroke;
    r = Cyl.CompressionRatio;  % Compression ratio
    
    % Calculate displacement volume Vd
    Vd = (pi / 4) * Bore^2 * Stroke;
    
    % Clearance volume Vc 
    Vc = Vd / (r - 1);
    
    % Volumes at key points
    V1 = Vc + Vd;  % Volume at Bottom Dead Center (BDC)
    V2 = Vc;       % Volume at Top Dead Center (TDC)
    
    % Process 1-2: Isentropic compression (from V1 to V2)
    T2 = T1 * (V1 / V2)^(k - 1);
    P2 = P1 * (V1 / V2)^k;
    
    % Volume at point 3 (after heat addition)
    V3 = rc * V2;
    
    % Process 2-3: Constant-pressure heat addition (from V2 to V3)
    P3 = P2;  % Pressure remains constant
    T3 = T2 * (V3 / V2);  % T3 = T2 * (V3 / V2)
    
    % Process 3-4: Isentropic expansion (from V3 to V4, where V4 = V1)
    V4 = V1;
    %T4_check = T3 * (V4 / V3)^(k - 1);  % Should match T4
    P4 = P3 * (V3 / V4)^k;
    
    % Verify T4
    %if abs(T4 - T4_check) > 1e-3
      %  warning('Calculated T4 does not match the given T4.');
    %end
    
    % Generate arrays for each process
    % Process 1-2: Isentropic compression
    V_12 = linspace(V1, V2, numPoints);
    P_12 = P1 * (V1 ./ V_12).^k;
    
    % Process 2-3: Constant-pressure heat addition
    V_23 = linspace(V2, V3, numPoints);
    P_23 = P2 * ones(size(V_23));
    
    % Process 3-4: Isentropic expansion
    V_34 = linspace(V3, V4, numPoints);
    P_34 = P3 * (V3 ./ V_34).^k;
    
    % Process 4-1: Constant-volume heat rejection (from P4 to P1)
    V_41 = V4 * ones(1, numPoints);
    P_41 = linspace(P4, P1, numPoints);
    
    % Concatenate the arrays to complete the cycle
    V_cycle = [V_12, V_23, V_34, V_41];
    P_cycle = [P_12, P_23, P_34, P_41];
end
