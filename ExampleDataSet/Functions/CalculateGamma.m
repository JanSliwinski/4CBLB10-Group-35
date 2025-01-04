function gamma = CalculateGamma(SpS,RPM,volume,pressure,O2percent,CO2percent,massFuel)
    %% Initialisation and checking if there is information on exhaust.
    Rair = 287; Runiv = 8.314; O2volfrac = O2percent/100; CO2volfrac = CO2percent/100;
    %% Exhuast composition data
    N2volfrac = 0.78;          % Nitrogen fraction in air
    O2volfrac_air = 0.21;      % Oxygen fraction in air
    H2Ovolfrac = O2volfrac_air - O2volfrac - CO2volfrac; % Calculate water vapor fraciton
    % Grouping together
    moleFractions = [N2volfrac, O2volfrac, CO2volfrac, H2Ovolfrac, 0]; % no diesel
    
    % Getting AFR
    AFR = compute_AFR(CO2volfrac,0, O2volfrac/100,N2volfrac);
    massFuel = mean(massFuel)/1000 / (RPM/60/2);
    massAir = massFuel * AFR;
    
    %% Compute specific gas constants for inlet and outlet
    % Get molar masses of species in kg/mol
    molarMass = [SpS.Mass];
    % Find mass fractions
    massFractions = moleFractions .* molarMass ./ (moleFractions * molarMass');
    % Find specific gas constants for all constituents
    Rs = Runiv./molarMass;
    % Get specific gas constats
    Rin = (Rair * AFR + Rs(5)) / (AFR+1);
    Rout = massFractions * Rs';
    
    %% Getting Temperature
    % Find maximum pressure and combustion point
    [~,combustion] = max(pressure);
    
    temperature = zeros(size(pressure));
    gamma = zeros(size(pressure));
    % %Temperature and gamma calculation loop
    for dummy = 1:size(gamma)
        % Pre-combustion calculations
        if dummy <= combustion
            % Calculate temperature using ideal gas law
            temperature(dummy) = pressure(dummy)*1e5 * volume(dummy) / (Rin * massAir);

            % Calculate specific heat ratio (gamma)
            gamma(dummy) = compute_cp(temperature(dummy),SpS,massFractions) ...
                / (compute_cp(temperature(dummy),SpS,massFractions) - Rin);

            % Post-combustion calculations
        else
            % Calculate temperature using exhaust gas properties
            temperature(dummy) = pressure(dummy)*1e5 * volume(dummy) / (Rout * massAir);

            % Calculate specific heat ratio (gamma)
            gamma(dummy) = compute_cp(temperature(dummy),SpS,massFractions) ...
                / (compute_cp(temperature(dummy),SpS,massFractions) - Rout);
        end
    end
    %% Vectorized temperature/gamma calculation (no compute_cp call)
% 
% nPoints    = length(pressure);
% idxPreComb = 1:combustion;
% idxPostComb= combustion+1 : nPoints;
% 
% % 1) Build a full temperature array (pre- or post-comb)
% temperature = zeros(nPoints,1);
% 
% % Pre-combustion: T = p*V/(Rin*massAir)
% temperature(idxPreComb) = ...
%     pressure(idxPreComb)*1e5 .* volume(idxPreComb) ./ (Rin  * massAir);
% 
% % Post-combustion: T = p*V/(Rout*massAir)
% temperature(idxPostComb) = ...
%     pressure(idxPostComb)*1e5 .* volume(idxPostComb) ./ (Rout * massAir);
% 
% % 2) Compute cp for *every point* and *every species* at once
% %    We'll store all species’ cp in one matrix. Each column = one species
% nSp      = numel(SpS);
% cpMatrix = zeros(nPoints, nSp);
% for i = 1:nSp
%     [cpVals, ~]       = CpNasa(temperature, SpS(i)); 
%     cpMatrix(:, i)    = cpVals;  
% end
% 
% % 3) Weighted sum of columns by massFractions => mixture cp at each point
% cpMix = cpMatrix * massFractions(:);  % size(cpMix) = [nPoints, 1]
% 
% % 4) Build an R array for entire domain (Rin vs. Rout)
% Rarray = zeros(nPoints,1);
% Rarray(idxPreComb)  = Rin;
% Rarray(idxPostComb) = Rout;
% 
% % 5) Finally, compute gamma = cp / (cp - R) for each point
% gamma = cpMix ./ (cpMix - Rarray);
end