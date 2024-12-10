function [Q_combustion] = MASSHeatOfCombustion(p_filtered, V_all, Ca, ValveEvents, gamma)
    % HeatOfCombustion calculates the heat of combustion for a given fuel.
    %
    % Inputs:
    %   fuel_name  - Name of the fuel (unused in current calculations)
    %   P          - Pressure vector (same length as Ca)
    %   V          - Volume vector (same length as Ca)
    %   Ca         - Crank angle vector
    %   FuelFR     - Fuel flow rate per unit crank angle - WHAT THE FLIP?
    %   ValveEvents- Structure containing CaIVC and CaEVO (start and end combustion angles)
    %   T_int      - Temp before combustion starts
    %
    % Output:
    %   Q_combustion - Calculated heat of combustion

    % Initialize Q_combustion
    Q_combustion = 0*size(Ca);
     
    % Ensure that P, V, Ca, and T_int are column vectors for consistency
    T_intcol = 295.15 * ones(3600,100);
    p_filtered = p_filtered(:);
    V_all = V_all(:);
    Ca = Ca(:);
    T_intcol = T_intcol(:);
    R = 8.314;
    % Validate input lengths
    if ~(length(p_filtered) == length(V_all) && length(V_all) == length(Ca)) == length(T_intcol)
        error('Inputs P, V, Ca and T_int must have the same length.');
    end
    StartCom = find(Ca == ValveEvents.CaIVC,1); % find CA where combustion starts
    EndCom = find(Ca == ValveEvents.CaEVO,1); % find CA where combustion ends

    %calculate amount of air in piston
    massair = (p_filtered(StartCom)*V_all(StartCom))/(R*T_intcol); %calculating mass of air in the cylinder before combustion by an ideal gas law
    
    % Compute gradients
    dVdCa = gradient(V_all, Ca);   % dV/dCa
    dPdCa = gradient(p_filtered, Ca);  %dP/dCa
    
    % Compute gamma using CpNasa
    % gammafuel = CpNasa(T_intcol, SpS)/CvNasa(T_intcol, SpS);    % Assuming CpNasa returns gamma for each T
    % gammaair = 1.4; % assming gamma_air = 1.4

    % Calculate Q_combustion using vectorized operations for efficiency
    for i = 2:length(Ca)
        if StartCom < Ca(i) && Ca(i) < EndCom
           % % OLD: massfuel = massfuel + FuelFR; % mass of fuel in engine?
           % massfuel = sum(mfr_fuel(i), 1) * Delta_t_perCA; % mass of fuel in engine
           %  totmass = massfuel+massair;   % total mass in engine - wouldn't this be the same as mass of gas reactants got from the ideal gas law?
           %  gamma = gammaair*massair/totmass + gammafuel*massfuel/totmass;
            %assume that theres no heat loss
            Q_combustion = Q_combustion(i-1) + gamma/(gamma-1) .* p_filtered(i)*dVdCa(i) + 1/(gamma-1) .* V_all(i) .* dPdCa(i); %apparent heat release /CA
        end
    end