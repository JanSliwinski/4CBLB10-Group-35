function KPITable = GenerateKPITable(KPIdataFiles, table2_experiment1, LHV, avg_m_fuelpercycle, RPM, AFR_stoich, x, MW_Fuel, Cyl)
%GENERATEKPITABLE Loads multiple datasets, calculates KPIs, and generates a summary table.
%
% Inputs:
%   KPIdataFiles       : Cell array with paths, fuel types, and crank angles.
%   table2_experiment1 : Table with emissions data (loaded from Excel in main script).
%   LHV                : Lower Heating Value of the fuel (J/kg).
%   avg_m_fuelpercycle : Mass of the fuel per cycle (g)
%   RPM                : Engine RPM.
%   AFR_stoich         : Stoichiometric Air-Fuel Ratio.
%   x                  : Coefficient of carbon in fuel (CxHy).
%   MW_Fuel            : Molecular weight of the fuel (g/mol).
%   Cyl                : Struct containing cylinder data
%
% Output:
%   KPITable           : Table containing calculated KPIs for each dataset.

    %% Initialize KPI Table
    n_rows = 7;
    KPITable = table(...
    strings(n_rows, 1), ...          % FuelType as string array
    zeros(n_rows, 1), ...            % CrankAngle
    zeros(n_rows, 1), ...            % Work
    zeros(n_rows, 1), ...            % Efficiency
    zeros(n_rows, 1), ...            % BSFC
    zeros(n_rows, 1), ...            % bsCO2
    zeros(n_rows, 1), ...            % bsNOx
    'VariableNames', {'FuelType', 'CrankAngle','Work', 'Efficiency', 'BSFC', 'bsCO2', 'bsNOx'});

    %% Loop Through Each Data File to Calculate KPIs
    for i = 1:size(KPIdataFiles, 1)
        % Extract metadata
        data_file_name = KPIdataFiles{i, 1};
        fuel_type = KPIdataFiles{i, 2};
        crank_angle = KPIdataFiles{i, 3};

        % Load Data
        data_in = table2array(readtable(data_file_name));

        % Reshape Data
        resolution = 0.2;  % Degrees crank angle resolution
        n_datapoints_per_cycle = 720 / resolution;
        n_cycles = size(data_in, 1) / n_datapoints_per_cycle;

        ca = reshape(data_in(:, 1), [], n_cycles);          % Crank angle in degrees
        p = reshape(data_in(:, 2), [], n_cycles) * 1e5;     % Pressure in Pa
        mfr_fuel = 0.16;    % Fuel mass flow rate (kg/s)

        % Apply Savitzky-Golay filter to pressure data
        polynomial_order = 3; % Adjust based on noise level
        frame_length = 21;    % Must be odd
        
        % Initialize the filtered pressure matrix
        p_filtered = zeros(size(p));
        
        % Apply the filter to each cycle
        for j = 1:n_cycles
            p_filtered(:, j) = SGFilter(p(:, j), polynomial_order, frame_length, 0);
        end
        
        % Replace raw pressure data with filtered data for further calculations
        p = p_filtered;

        % Calculate Work and Power
        v_all = zeros(size(ca));
        for j = 1:n_cycles
            v_all(:, j) = CylinderVolume(ca(:, j), Cyl);
        end

        v_avg = mean(v_all, 2);
        p_avg = mean(p, 2);
        W = trapz(v_avg, p_avg);             % Work done (J)
        P = W * (RPM / (2 * 60));            % Power output (W)

        % Calculate Air Mass Flow Rate
        o2_percent_load = table2_experiment1{:, 'O2_percent'};
        mfr_air = CalculateMassFlowAir(o2_percent_load(i), mfr_fuel, AFR_stoich);

        % Get NOx emissions for this crank angle
        target_crank_angle = crank_angle;
        row_idx = find(table2_experiment1{:, 'CrankAngle'} == target_crank_angle);

        if ~isempty(row_idx)
            nox_ppm = table2_experiment1{row_idx, 'NOx_ppm'};
        else
            error(['NOx data for CrankAngle ', num2str(target_crank_angle), ' not found.']);
        end

        % Calculate KPIs
        KPIs = CalculateKPIs(W, mean(mfr_fuel, 'all'), avg_m_fuelpercycle, LHV, P, x, mean(mfr_air, 'all'), nox_ppm, MW_Fuel);

        % Populate the i-th row of KPITable
        KPITable(i, :) = {fuel_type, crank_angle, W, KPIs.Efficiency, KPIs.BSFC, KPIs.bsCO2, KPIs.bsNOx};
    end

end
