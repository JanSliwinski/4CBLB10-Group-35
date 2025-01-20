function [Q_combustion_aROHR] = Calc_Q_aROHR(p_filtered, V_all, Ca, ValveEvents, gamma)
    % HeatOfCombustion calculates the heat of combustion for a given fuel.
    %
    % Inputs:
    %   fuel_name  - Name of the fuel (unused in current calculations)
    %   P          - Pressure vector (same length as Ca)
    %   V          - Volume vector (same length as Ca)
    %   Ca         - Crank angle vector
    %   ValveEvents- Structure containing CaIVC and CaEVO (start and end combustion angles)
    %   T_int      - Temp before combustion starts
    %
    % Output:
    %   Q_combustion - Calculated heat of combustion
     
    % Ensure that P, V, Ca, and T_int are column vectors for consistency
    T_intcol = 295.15 * ones(size(Ca));
    p_filtered = p_filtered(:);
    V_all = V_all(:);
    Ca = Ca(:);
    T_intcol = T_intcol(:);
    R = 8.314;
    % Validate input lengths
    if ~(length(p_filtered) == length(V_all) && length(V_all) == length(Ca)) == length(T_intcol)
        error('Inputs P, V, Ca and T_int must have the same length.');
    end
    % Initialize Q_combustion
    Q_combustion = zeros(size(Ca));

    StartCom = find(Ca == ValveEvents.CaIVC,1); % find CA where combustion starts
    EndCom = find(Ca == ValveEvents.CaEVO,1); % find CA where combustion ends

    %calculate amount of air in piston
    massair = (p_filtered(StartCom)*V_all(StartCom))/(R*T_intcol); %calculating mass of air in the cylinder before combustion by an ideal gas law
    
    % Compute gradients
    dVdCa = gradient(V_all, Ca);   % dV/dCa
    dPdCa = gradient(p_filtered, Ca);  %dP/dCa

    %verify that the indices of start and end exist
    disp(StartCom);
    disp(EndCom);
    % Calculate Q_combustion using vectorized operations for efficiency
    for i = 2:length(Ca)
    if i >= StartCom && i <= EndCom      
        % Update Q_combustion
        Update = gamma/(gamma-1) * p_filtered(i) * dVdCa(i) + 1/(gamma-1) * V_all(i) * dPdCa(i);
        Q_combustion(i) = Q_combustion(i-1) + Update;
    end
    Q_combustion_isolated = Q_combustion(StartCom:EndCom);
    Q_combustion_aROHR = mean(Q_combustion_isolated);
    % 
    % figure;
    % plot(Ca(StartCom:EndCom), Q_combustion_isolated, 'LineWidth', 1.5);  % Ca(:,1) is the crank angle array
    % xlabel('Crank Angle (°)');
    % ylabel('Energy of combustion [J]');
    % title('Energy of combustion (isolated)');
    % %xlim([-45,135]);
    % grid on;
    % hold on;
    end
   