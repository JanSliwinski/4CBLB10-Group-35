function gamma = CalculateGamma(T,expNo)    
    %% Initialisation and checking if there is information on exhaust.
    CA = (0:0.2:720-0.2)'; RPM=1500; Runiv = 8.314; mm=0.001; Rair = 287;

    if isstruct(T.AdditionalData{expNo})
        %% Load Nasa database
        addpath("Nasa\", "Functions\");
        TdataBase = fullfile('Nasa', 'NasaThermalDatabase'); 
        load(TdataBase);
        
        %% Find species
        iSp = myfind({Sp.Name}, {'N2', 'O2', 'CO2','H2O','Diesel'}); % Find indexes of these species
        SpS = Sp(iSp); % Create subset of species based on found indexes
        
        %% Exhuast composition data
        N2volfrac = 0.78;          % Nitrogen fraction in air
        O2volfrac_air = 0.21;      % Oxygen fraction in air
        O2volfrac = T.AdditionalData{expNo}.O2 / 100; % Oxygen fraction
        CO2volfrac = T.AdditionalData{expNo}.CO2 / 100; % Carbon dioxide fracitons
        H2Ovolfrac = O2volfrac_air - O2volfrac - CO2volfrac; % Calculate water vapor fraciton
        % Grouping together
        moleFractions = [N2volfrac, O2volfrac, CO2volfrac, H2Ovolfrac, 0]; % no diesel
    
        % Getting AFR
        AFR = compute_AFR(CO2volfrac,0, O2volfrac,N2volfrac);
        massFuel = mean(T.AverageCycleData{expNo}.AvgMassFlow)/1000 / (RPM/60/2);
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
        
        %% getting air mass per cycle and cylinder volume
        % Engine Geometry Parameters
        Cyl.Bore = 104 * mm;               % Cylinder bore
        Cyl.Stroke = 85 * mm;              % Cylinder stroke
        Cyl.CompressionRatio = 21.5;       % Compression ratio
        Cyl.ConRod = 136.5 * mm;           % Connecting rod length
        Cyl.TDCangle = 180;                % Top Dead Center angle  
        % Calculate cylinder volume using CylinderVolume function
        volume = CylinderVolume(CA,Cyl);
    
        %% Getting Temperature
        % Find maximum pressure and combustion point
        pressure = T.AverageCycleData{expNo}.AvgPressure;
        [~,combustion] = max(pressure);
        
        temperature = zeros(size(pressure));
        gamma = zeros(size(pressure));
        % Temperature and gamma calculation loop
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
        
    else
        gamma = 1.32; % Default value when empty (NaN)
    end
end