function gamma = CalculateGamma(SpS, RPM, volume, pressure, O2percent, CO2percent, massFuel)
    %% CalculateGamma: Computes the ratio of specific heats (gamma) at each point
    % in the crank-angle domain. Pre- and post-combustion are handled separately,
    % and a mixture of gases is considered using NASA polynomial data for cp.

    %% Initial setup
    Rair  = 287;       % Specific gas constant for air [J/(kg·K)]
    Runiv = 8.314;     % Universal gas constant [J/(mol·K)]
    % Convert the O2 and CO2 input (in percent) to fractional volumes
    O2volfrac  = O2percent  / 100; 
    CO2volfrac = CO2percent / 100;

    %% Exhaust composition data
    N2volfrac      = 0.78;          % Approx. nitrogen fraction in air
    O2volfrac_air  = 0.21;          % Oxygen fraction in air
    % Water vapor fraction: difference between air O2 fraction and
    % the measured O2 + CO2 in exhaust
    H2Ovolfrac     = O2volfrac_air - O2volfrac - CO2volfrac;
    % Combine all relevant species' volume fractions
    moleFractions  = [N2volfrac, O2volfrac, CO2volfrac, H2Ovolfrac, 0]; % The last '0' is for Diesel

    %% Calculate Air-Fuel Ratio (AFR)
    % compute_AFR() is presumably a custom function that uses O2, CO2, etc.
    AFR = compute_AFR(CO2volfrac, 0, O2volfrac / 100, N2volfrac);

    % Convert mass flow from mg/s to kg/cycle:
    % - mean(...) picks the average injection flow rate
    % - dividing by 1000 converts mg to g
    % - dividing by (RPM/60/2) converts from grams/s to grams/cycle 
    %   for a 4-stroke cycle at the given RPM
    massFuel = mean(massFuel) / 1000 / (RPM / 60 / 2);

    % Estimate mass of air based on the mass of fuel and AFR
    massAir = massFuel * AFR;

    %% Compute gas constants (R) for mixture
    % Extract species' molar masses (kg/mol) from the NASA database entries (SpS)
    molarMass    = [SpS.Mass];
    % Convert volume fractions to mass fractions:
    massFractions = moleFractions .* molarMass ./ (moleFractions * molarMass');

    % Compute each species' individual gas constant = Runiv / (molar mass)
    Rs = Runiv ./ molarMass;

    % Inlet mixture gas constant Rin: a combination of air + diesel
    % Weighted by the AFR
    Rin  = (Rair * AFR + Rs(5)) / (AFR + 1); 
    % Exhaust mixture gas constant Rout: a mass-weighted combination 
    % of the species present post-combustion
    Rout = massFractions * Rs';

    %% Determine combustion point (peak pressure)
    % The index at which pressure is max is considered the "combustion" boundary
    [~, combustion] = max(pressure);

    %% Prepare to compute temperature and gamma across the crank angle range
    nPoints    = length(pressure);  % Total number of points in the domain
    idxPreComb = 1 : combustion;    % Indices before combustion peak
    idxPostComb= combustion+1 : nPoints; % Indices after combustion peak

    % Preallocate arrays
    temperature = zeros(nPoints, 1);

    %% Compute Temperature
    % Pre-combustion: T = p * V / (Rin * massAir), using inlet mixture R
    temperature(idxPreComb) = ...
        pressure(idxPreComb) * 1e5 .* volume(idxPreComb) ./ (Rin  * massAir);

    % Post-combustion: T = p * V / (Rout * massAir), using exhaust mixture R
    temperature(idxPostComb) = ...
        pressure(idxPostComb)* 1e5 .* volume(idxPostComb) ./ (Rout * massAir);

    %% Compute cp for each point
    % Build a matrix to hold cp values for each species at every crank-angle point
    nSp      = numel(SpS);         % Number of species
    cpMatrix = zeros(nPoints, nSp);% Will store cp for each species across all points

    % Loop over species, compute cp with NASA polynomials at each temperature
    for i = 1 : nSp
        [cpVals, ~]   = CpNasa(temperature, SpS(i)); 
        cpMatrix(:, i)= cpVals;  
    end

    % Weighted sum of species cp => mixture cp
    % Each row of cpMatrix is a set of species cp at a certain crank angle.
    % Multiplying by massFractions(:) sums them up by their mass fraction.
    cpMix = cpMatrix * massFractions(:);

    %% Build an R array for entire domain (Rin vs. Rout)
    Rarray            = zeros(nPoints, 1);
    Rarray(idxPreComb)= Rin;
    Rarray(idxPostComb) = Rout;

    %% Compute gamma = cp / (cp - R) at every point
    gamma = cpMix ./ (cpMix - Rarray);

end
