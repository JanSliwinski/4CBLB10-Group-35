function KPIs = CalculateKPIs(W, mfr_fuel, avg_m_fuelpercycle, LHV, P, x, mfr_air, NOx_ppm, MW_Fuel)
% Inputs:
%   x                  : Coefficient of carbon in fuel (CxHy)
%   NOx_ppm_vector     : Vector with NOx part per million reading for each CA (ppm)
%   W                  : Work done Vector with each work done for each CA (J)
%   avg_m_fuelpercycle : Mass of the fuel per cycle (g)
%   mfr_fuel           : Fuel mass flow rate(kg/s) (Should be the same as same load at each CA)
%   LHV                : Lower heating value of the fuel (J/kg)
%   P_vector           : Power output for each  CA (W)
%   mfr_air            : Air mass flow rate (kg/s)

% Outputs:
%   KPIs     : Struct containing calculated KPIs

MW_CO2 = 44.01;


%% eg
% Input data (replace these with data from the Excel file)
CO2_frac = 4.7 / 100;  % Mole fraction of CO2 (from percentage)
O2_frac = 14.35 / 100; % Mole fraction of O2 (from percentage)
N2_frac = 1 - (CO2_frac + O2_frac); % Assume the remaining is N2
NOx_frac = 377 / 1e6;  % Mole fraction of NOx (from ppm)

% Molar masses (g/mol)
M_CO2 = 44;
M_O2 = 32;
M_N2 = 28;
M_NOx = 38; % Average molar mass of NOx
M_exhaust = (CO2_frac * M_CO2) + (O2_frac * M_O2) + (N2_frac * M_N2) + (NOx_frac * M_NOx);




% Efficiency
KPIs.Efficiency = P / (0.16 * LHV);

% Brake-Specific Fuel Consumption (BSFC)
KPIs.BSFC = (mfr_fuel  * 3600) ./ (P / 1000); % g/kWh

% Brake-Specific CO2 Emissions
KPIs.bsCO2 = ((x * MW_CO2) / MW_Fuel) * KPIs.BSFC;

% Brake-Specific NOx Emissions
mfr_NOx = (NOx_ppm / 1e6)*(38/M_exhaust) * (mfr_air + 0.16); % Assuming exhaust mass flow is total
KPIs.bsNOx = mfr_NOx / (P/1000);

end
