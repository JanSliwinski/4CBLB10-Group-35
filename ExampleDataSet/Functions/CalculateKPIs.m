function KPIs = CalculateKPIs(true_mfr_fuel, LHV, P, x, mfr_air, NOx_ppm, MW_Fuel)
    % Inputs:
    %   true_mfr_fuel      : Fuel mass flow rate (kg/s)
    %   LHV                : Lower heating value of the fuel (J/kg)
    %   P                  : Power output (W)
    %   x                  : Number of carbon atoms in the fuel (CxHy)
    %   mfr_air            : Air mass flow rate (kg/s)
    %   NOx_ppm            : NOx concentration in ppm (by volume)
    %   MW_Fuel            : Molar mass of the fuel (g/mol)

    % Outputs:
    %   KPIs : A struct containing:
    %          - Efficiency
    %          - BSFC (g/kWh)
    %          - bsCO2 (g/kWh)
    %          - bsNOx (g/kWh)

    % Define constants
    MW_CO2 = 44.01;   % g/mol CO2
    MW_NOx = 46;      % g/mol NO2 (used as reference for NOx)
    
    % Assumed exhaust gas composition for calculating M_exhaust
    % Ideally, these should come from measurements or given data.
    CO2_frac = 4.7 / 100;    % Mole fraction of CO2
    O2_frac  = 14.35 / 100;  % Mole fraction of O2
    N2_frac  = 1 - (CO2_frac + O2_frac); 

    M_CO2 = 44;  % g/mol
    M_O2  = 32;  % g/mol
    M_N2  = 28;  % g/mol

    % Convert NOx from ppm to mole fraction
    NOx_frac = NOx_ppm / 1e6;

    % Mean molecular weight of exhaust (g/mol)
    M_exhaust = (CO2_frac * M_CO2) + (O2_frac * M_O2) + (N2_frac * M_N2) + (NOx_frac * MW_NOx);

    % Efficiency
    % Check reference factor (0.16) as needed
    KPIs.Efficiency = P / (true_mfr_fuel * LHV);

    % Brake-Specific Fuel Consumption (BSFC)
    % (true_mfr_fuel [g/s]*3600 [s/h])/(P [W]/1000 [kW]) = g/kWh
    KPIs.BSFC = ((true_mfr_fuel * 3600) / (P / 1000));

    % Brake-Specific CO2 Emissions
    % bsCO2 = ((x * MW_CO2)/MW_Fuel)*BSFC
    KPIs.bsCO2 = ((x * MW_CO2) / MW_Fuel) * KPIs.BSFC;

    % Brake-Specific NOx Emissions
    % Convert NOx ppm to mass fraction
    NOx_mass_fraction = NOx_frac * (MW_NOx / M_exhaust);

    % Total exhaust mass flow (approx. sum of air and fuel)
    mfr_exhaust = mfr_air + true_mfr_fuel;

    % NOx mass flow (kg/s)
    mfr_NOx = NOx_mass_fraction * mfr_exhaust;

    % Convert NOx mass flow to g/kWh
    % (mfr_NOx [g/s]*3600 [s/h])/(P [W]/1000 [kW]) = g/kWh
    KPIs.bsNOx = ((mfr_NOx * 3600) / (P / 1000));
end
