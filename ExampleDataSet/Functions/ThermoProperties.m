 function [Delta_H_all, Delta_U_all, Delta_S_all] = ThermoProperties(T_int, T_exh, SpS, Ncycles, Ca, intake_species, exhaust_species, Y_int, Y_exh)
    % Inputs:
    % T_int           : Matrix of intake temperatures (rows = crank angles, cols = cycles)
    % T_exh           : Matrix of exhaust temperatures (rows = crank angles, cols = cycles)
    % SpS             : NASA polynomial coefficients for species
    % Ncycles         : Number of engine cycles
    % Ca              : Crank angle data matrix (used for sizing outputs)
    % intake_species  : Indexes of species present in the intake 
    % exhaust_species : Indexes of species present in the exhaust 
    % Y_int           : Mass or mole fractions of intake species 
    % Y_exh           : Mass or mole fractions of exhaust species
    % 
    % Outputs:
    % Delta_H_all     : Matrix of enthalpy changes for all cycles
    % Delta_U_all     : Matrix of internal energy changes for all cycles
    % Delta_S_all     : Matrix of entropy changes for all cycles

    % Initialize output matrices
    Delta_H_all = zeros(size(Ca));
    Delta_U_all = zeros(size(Ca));
    Delta_S_all = zeros(size(Ca));

    % Loop through each cycle to compute properties
    for i = 1:Ncycles
        % Initializes empty vectors for intake
        H_int = zeros(size(T_int(:, i)));
        U_int = zeros(size(T_int(:, i)));
        S_int = zeros(size(T_int(:, i)));
        
        % Initializes empty vectors for exhaust
        H_exh = zeros(size(T_exh(:, i)));
        U_exh = zeros(size(T_exh(:, i)));
        S_exh = zeros(size(T_exh(:, i)));
        
        % Sum contributions from all intake species
        for idx = intake_species
            H_int = H_int + Y_int * HNasa(T_int(:, i), SpS(idx));  % Weighted Enthalpy
            U_int = U_int + Y_int * UNasa(T_int(:, i), SpS(idx));  % Weighted Internal energy
            S_int = S_int + Y_int * SNasa(T_int(:, i), SpS(idx));  % Weighted Entropy
        end
        
        % Sum contributions from all exhaust species
        for idx = exhaust_species
            H_exh = H_exh + Y_exh(idx) * HNasa(T_exh(:, i), SpS(idx));  % Weighted Enthalpy
            U_exh = U_exh + Y_exh(idx) * UNasa(T_exh(:, i), SpS(idx));  % Weighted Internal energy
            S_exh = S_exh + Y_exh(idx) * SNasa(T_exh(:, i), SpS(idx));  % Weighted Entropy
        end

        % Compute changes
        Delta_H_all(:, i) = H_exh - H_int;
        Delta_U_all(:, i) = U_exh - U_int;
        Delta_S_all(:, i) = S_exh - S_int;
    end
end
