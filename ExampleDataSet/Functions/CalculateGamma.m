function gamma = CalculateGamma(SpS, volume, pressure, O2percent, CO2percent, mfr_fuel, AFR, RPM)
    try
        %% Input validation
        validateInput = @(x, name) assert(~isempty(x) && isnumeric(x) && all(isfinite(x)), ...
            'Error: %s must be non-empty, numeric, and finite', name);
        
        % Validate basic numeric inputs
        validateInput(volume, 'volume');
        validateInput(pressure, 'pressure');
        validateInput(O2percent, 'O2percent');
        validateInput(CO2percent, 'CO2percent');
        validateInput(mfr_fuel, 'true_mfr_fuel');
        validateInput(AFR, 'AFR_stoich');
        validateInput(RPM, 'RPM');
        
        % Validate SpS structure
        assert(isstruct(SpS) && isfield(SpS, 'Mass'), 'Error: SpS must be a structure with Mass field');
        validateInput([SpS.Mass], 'SpS.Mass');
        
        % Validate array lengths
        assert(length(volume) == length(pressure), ...
            'Error: volume and pressure arrays must have the same length');
        
        % Validate ranges
        assert(all(O2percent >= 0 & O2percent <= 100), 'Error: O2percent must be between 0 and 100');
        assert(all(CO2percent >= 0 & CO2percent <= 100), 'Error: CO2percent must be between 0 and 100');
        assert(all(RPM > 0), 'Error: RPM must be positive');
        assert(all(AFR > 0), 'Error: AFR_stoich must be positive');
        assert(all(mfr_fuel > 0), 'Error: true_mfr_fuel must be positive');
        
        %% Constants
        Rair = 287;   % Specific gas constant for air [J/(kg·K)]
        Runiv = 8.314; % Universal gas constant [J/(mol·K)]
        
        %% Convert percentages to fractions with safety checks
        O2volfrac = O2percent / 100;
        CO2volfrac = CO2percent / 100;
        
        %% Exhaust composition calculations
        N2volfrac = 0.78;
        O2volfrac_air = 0.21;
        
        % Check for physically realistic water vapor fraction
        H2Ovolfrac = O2volfrac_air - O2volfrac - CO2volfrac;
        assert(H2Ovolfrac >= 0, ['Error: Calculated H2O volume fraction is negative. ', ...
            'Check O2 and CO2 percentages.']);
        
        moleFractions = [N2volfrac, O2volfrac, CO2volfrac, H2Ovolfrac, 0];
        assert(abs(sum(moleFractions) - 1) < 0.1, ...
            'Error: Mole fractions sum deviates significantly from 1.0');
        
        %% Mass calculations with unit conversion safety
        cyclesPerSecond = RPM / 60 / 2; % For 4-stroke engine
        assert(cyclesPerSecond > 0, 'Error: Invalid RPM leads to zero or negative cycles per second');
        
        % Convert fuel mass flow from mg/s to kg/cycle
        massAir = (mfr_fuel/1000) / cyclesPerSecond * AFR;
        assert(massAir > 0, 'Error: Calculated air mass is not positive');
        
        %% Gas constant calculations
        molarMass = [SpS.Mass];
        assert(all(molarMass > 0), 'Error: All molar masses must be positive');
        
        % Mass fraction calculations
        denominator = moleFractions * molarMass';
        assert(denominator > 0, 'Error: Invalid denominator in mass fraction calculation');
        massFractions = moleFractions .* molarMass ./ denominator;
        
        % Individual gas constants
        Rs = Runiv ./ molarMass;
        
        % Mixture gas constants
        Rin = (Rair * AFR + Rs(5)) / (AFR + 1);
        Rout = massFractions * Rs';
        
        assert(Rin > 0 && Rout > 0, 'Error: Invalid mixture gas constants calculated');
        
        %% Combustion point determination
        [~, combustion] = max(pressure);
        assert(combustion > 1 && combustion < length(pressure), ...
            'Error: Combustion point at array boundary - check pressure data');
        
        %% Temperature calculations
        nPoints = length(pressure);
        idxPreComb = 1:combustion;
        idxPostComb = (combustion+1):nPoints;
        
        temperature = zeros(nPoints, 1);
        
        % Pre-combustion temperatures
        temperature(idxPreComb) = pressure(idxPreComb) * 1e5 .* volume(idxPreComb) / (Rin * massAir);
        
        % Post-combustion temperatures
        temperature(idxPostComb) = pressure(idxPostComb) * 1e5 .* volume(idxPostComb) / (Rout * massAir);
        
        % Validate temperatures
        assert(all(temperature > 0), 'Error: Invalid negative or zero temperatures calculated');
        assert(all(isfinite(temperature)), 'Error: Non-finite temperatures calculated');
        
        %% Specific heat calculations
        nSp = numel(SpS);
        cpMatrix = zeros(nPoints, nSp);
        
        % Calculate cp for each species
        for i = 1:nSp
            try
                [cpVals, ~] = CpNasa(temperature, SpS(i));
                assert(all(isfinite(cpVals)), 'Non-finite cp values calculated');
                cpMatrix(:, i) = cpVals;
            catch ME
                error('Error calculating cp for species %d: %s', i, ME.message);
            end
        end
        
        % Calculate mixture cp
        cpMix = cpMatrix * massFractions(:);
        assert(all(cpMix > 0), 'Error: Invalid negative or zero specific heats calculated');
        
        %% Final gamma calculation
        Rarray = zeros(nPoints, 1);
        Rarray(idxPreComb) = Rin;
        Rarray(idxPostComb) = Rout;
        
        gamma = cpMix ./ (cpMix - Rarray);
        
        % Validate final gamma values
        assert(all(gamma > 1) && all(gamma < 2), ...
            'Error: Gamma values outside physically reasonable range (1 < γ < 2)');
        assert(all(isfinite(gamma)), 'Error: Non-finite gamma values calculated');
        
    catch ME
        % Enhanced error reporting
        errorMsg = sprintf(['Error in CalculateGamma:\n', ...
            'Message: %s\n', ...
            'Function: %s\n', ...
            'Line: %d\n'], ...
            ME.message, ME.stack(1).name, ME.stack(1).line);
        
        % Log error to file
        logError(errorMsg);
        
        % Rethrow with additional context
        error('GammaCalculationError:Failed', '%s', errorMsg);
    end
end

function logError(errorMsg)
    % Create a log file with timestamp
    try
        logFile = 'gamma_calculation_errors.log';
        fid = fopen(logFile, 'a');
        if fid ~= -1
            fprintf(fid, '\n%s\n%s\n', datestr(now), errorMsg);
            fclose(fid);
        end
    catch
        % If logging fails, don't let it interrupt the error handling
        warning('Failed to log error to file');
    end
end