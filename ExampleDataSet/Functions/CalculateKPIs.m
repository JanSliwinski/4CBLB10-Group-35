function KPIs = CalculateKPIs(W, mfr_fuel, LHV, P, x, mfr_air, NOx_ppm, MW_Fuel)
% Inputs:
%   x        : Coefficient of carbon in fuel (CxHy)
%   NOx_ppm  : NOx part per million reading from experiment
%   W        : Work done (J)
%   mfr_fuel : Fuel mass flow rate (kg/s)
%   LHV      : Lower heating value of the fuel (J/kg)
%   P        : Power output (W)
%   mfr_air  : Air mass flow rate (kg/s)

% Outputs:
%   KPIs     : Struct containing calculated KPIs

MW_CO2 = 44.01;
% Efficiency
KPIs.Efficiency = W / (mfr_fuel * LHV);

% Brake-Specific Fuel Consumption (BSFC)
KPIs.BSFC = ((mfr_fuel * 1000) * 3600) / (P / 1000); % g/kWh

% Brake-Specific CO2 Emissions
KPIs.bsCO2 = ((x * MW_CO2) / MW_Fuel) * KPIs.BSFC;

% Brake-Specific NOx Emissions
mfr_NOx = (NOx_ppm / 10e6) * (mfr_air + mfr_fuel); % Assuming exhaust mass flow is total
KPIs.bsNOx = mfr_NOx / P;

end