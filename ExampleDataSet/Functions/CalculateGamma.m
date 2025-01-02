function gamma = CalculateGamma(T,expNo)
    Runiv = 8.314; Tref = 293.15; mm=0.001;
    
    %% Unpacking the table T
    CA = (0:0.2:720-0.2)';
    pressure = T.AverageCycleData{expNo}.AvgPressure;
    
    %% Load Nasa database
    addpath("Nasa\");
    TdataBase = fullfile('Nasa', 'NasaThermalDatabase'); 
    load(TdataBase);
    
    %% Find species
    iSp = myfind({Sp.Name}, {'O2', 'CO2', 'N2','H2O'}); % Find indexes of these species
    SpS = Sp(iSp); % Create subset of species based on found indexes
    
    %% Exhuast composition data
    N2_volfrac = 0.78;          % Nitrogen fraction in air
    O2_volfrac_air = 0.21;      % Oxygen fraction in air
    O2_volfrac = T.AdditionalData{expNo}.O2 / 100; % Oxygen fraction
    CO2_volfrac = T.AdditionalData{expNo}.CO2 / 100; % Carbon dioxide fracitons
    H2O_volfrac = O2_volfrac_air - O2_volfrac - CO2_volfrac; % Calculate water vapor fraciton by subtraction
    
    %% Compute specific gas constants for inlet and outlet
    % Get molar masses of species in kg/mol
    Mi = [SpS.Mass];
    % Find specific gas constants for all constituents
    R_co2 = Runiv/Mi(2); R_o2 = Runiv/Mi(1); R_n2 = Runiv/Mi(3); R_h20 = Runiv/Mi(4);
    
    % Find mass fractions
    totalMass = N2_volfrac*Mi(3)+O2_volfrac*Mi(1)+H2O_volfrac*Mi(4)+CO2_volfrac*Mi(2);
    N2_massfrac = N2_volfrac*Mi(3) / totalMass;
    O2_massfrac = O2_volfrac*Mi(1) / totalMass;
    CO2_massfrac = CO2_volfrac*Mi(2) / totalMass;
    H2O_massfrac = H2O_volfrac*Mi(4) / totalMass;
    massFractions = [O2_massfrac,CO2_massfrac,N2_massfrac,H2O_massfrac];
    
    % Get specific gas constats
    R_in = 287; % Ambient air
    R_out = N2_massfrac*R_n2+O2_massfrac*R_o2+CO2_massfrac*R_co2+H2O_massfrac*R_h20;
    
    %% getting air mass per cycle and cylinder volume
    % Engine Geometry Parameters
    Cyl.Bore = 104 * mm;               % Cylinder bore
    Cyl.Stroke = 85 * mm;              % Cylinder stroke
    Cyl.CompressionRatio = 21.5;       % Compression ratio
    Cyl.ConRod = 136.5 * mm;           % Connecting rod length
    Cyl.TDCangle = 180;                % Top Dead Center angle
    
    % Calculate cylinder volume using CylinderVolume function
    volume = CylinderVolume(CA,Cyl);
    
    [maxVolume,idx] = max(volume);
    massAir = pressure(idx)*1e5 * maxVolume / (R_in * Tref);
    
    %% Getting Temperature
    
    % Find maximum pressure and combustion point
    [~,combust] = max(pressure);
    
    temperature = zeros(size(pressure));
    gamma = zeros(size(pressure));
    % Temperature and gamma calculation loop
    for dummy = 1:size(gamma)
        % Pre-combustion calculations
        if dummy <= combust
            % Calculate temperature using ideal gas law
            temperature(dummy) = pressure(dummy)*1e5 * volume(dummy) / (R_in * massAir);
            
            % Calculate specific heat ratio (gamma)
            gamma(dummy) = compute_cp(temperature(dummy),SpS,massFractions) ...
                / (compute_cp(temperature(dummy),SpS,massFractions) - R_in);
        
        % Post-combustion calculations
        else
            % Calculate temperature using exhaust gas properties
            temperature(dummy) = pressure(dummy)*1e5 * volume(dummy) / (R_out * massAir);
            
            % Calculate specific heat ratio (gamma)
            gamma(dummy) = compute_cp(temperature(dummy),SpS,massFractions) ...
                / (compute_cp(temperature(dummy),SpS,massFractions) - R_out);
        end
    end
end