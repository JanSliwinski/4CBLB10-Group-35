function KPITable = GenerateKPITable(IDsforKPI, true_mfr_fuel, T, LHV, RPM, AFR_stoich, x, MW_fuel, Cyl, fuel_used)
%GENERATEKPITABLE Loads multiple datasets, calculates KPIs, and generates a summary table.
%
% Inputs:
%   KPIdataFiles       : Cell array with paths, fuel types, and crank angles.
%   true_mfr_fuel      : The calulated mass flow during the injection window(g/s)
%   table2_experiment1 : Table with emissions data (loaded from Excel in main script).
%   LHV                : Lower Heating Value of the fuel (J/kg).
%   RPM                : Engine RPM.
%   AFR_stoich         : Stoichiometric Air-Fuel Ratio.
%   x                  : Coefficient of carbon in fuel (CxHy).
%   M_fuel             : Molecular weight of the fuel (g/mol).
%   Cyl                : Struct containing cylinder data
%
% Output:
%   KPITable           : Table containing calculated KPIs for each dataset.

    %% Initialize KPI Table
    n_rows = height(IDsforKPI);
    KPITable = table(...
    strings(n_rows, 1), ...          % FuelType as string array
    zeros(n_rows, 1), ...            % CrankAngle
    zeros(n_rows, 1), ...            % Work
    zeros(n_rows, 1), ...            % Efficiency
    zeros(n_rows, 1), ...            % BSFC
    zeros(n_rows, 1), ...            % bsCO2
    zeros(n_rows, 1), ...            % bsNOx
    'VariableNames', {'FuelType', 'CrankAngle[°]','Work[J]', 'Efficiency[-]', 'BSFC[g/kWh]', 'bsCO2[g/kWh]', 'bsNOx[g/kWh]'});

   crank_angles = 14:20;


    %% Loop Through Each Data File to Calculate KPIs
    for i = 1:size(IDsforKPI)
        % Extract metadata
        data_file_name = IDsforKPI{i,1};
    
        % Extract fuel type (example logic, adjust as needed)
        fuel_type = fuel_used; % Example: HVO50
    
        % Load Data
        rowIndex = find(strcmp(T.UniqueID, data_file_name));   % Find the row index for the desired ID
        data_in = T.AverageCycleData{rowIndex};
        experimentDataCell = T.ExperimentData{rowIndex};    % Get the 1x1 cell of ExperimentData
        experimentData = experimentDataCell{1};         % Extract the 36000x4 matrix inside
        
        %hard coding some important parts:
        crank_angle = crank_angles(i);

        ca = experimentData(1:3600, 1);   % Crank angle in degrees
        p = data_in.AvgPressure * 1e5;     % Pressure in Pa
        mfr_fuel = true_mfr_fuel;    % Fuel mass flow rate (kg/s)

        % Apply Savitzky-Golay filter to pressure data
        polynomial_order = 3; % Adjust based on noise level
        frame_length = 21;    % Must be odd
        
        % Initialize the filtered pressure matrix
        p_filtered = zeros(size(p));
        
         % Apply the filter for one cycle
        p_filtered = SGFilter(p, polynomial_order, frame_length, 0);
        
        % Replace raw pressure data with filtered data for further calculations
        p = p_filtered;

        % Calculate Work and Power
        v_all = CylinderVolume(ca, Cyl);

        v_avg = mean(v_all, 2);
        p_avg = mean(p, 2);
        W = trapz(v_avg, p_avg);             % Work done (J)
        P = W * (RPM / (2 * 60));            % Power output (W)

        exhaustDatainT = T.AdditionalData{rowIndex};    % Extract the Exhaust data for the blend

        % Calculate Air Mass Flow Rate
        O2_percent = exhaustDatainT.O2;
        mfr_air = CalculateMassFlowAir(O2_percent, mfr_fuel, AFR_stoich);

        % Get emissions data
        nox_ppm = exhaustDatainT.NOx;
        % CO_percent = exhaustDatainT.CO;
        % HC_ppm = exhaustDatainT.HC;
        CO2_percent = exhaustDatainT.CO2;
        O2_percent = exhaustDatainT.O2;

        % Calculate KPIs
        KPIs = CalculateKPIs(true_mfr_fuel, LHV, P, x, mean(mfr_air, 'all'), nox_ppm, MW_fuel, CO2_percent, O2_percent);

        % Populate the i-th row of KPITable
        KPITable(i, :) = {fuel_type, crank_angle, W, KPIs.Efficiency, KPIs.BSFC, KPIs.bsCO2, KPIs.bsNOx};
    end
disp(mfr_air)
end
