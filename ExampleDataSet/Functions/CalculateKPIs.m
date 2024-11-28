function KPIs = CalculateKPIs(W, mfr_fuel, LHV, P, mfr_CO2, mfr_NOx)
% Inputs:
% W         : Work done (J)
% mfr_fuel  : Fuel mass flow rate (kg/s)
% LHV       : Lower heating value of the fuel (J/kg)
% P         : Power output (W)
% mfr_CO2   : CO2 mass flow rate (Kg/s)
% mfr_NOx   : NOx mass flow rate (Kg/s)
%
% Outputs:
% KPIs      : Struct containing calculated KPIs

    % Efficiency
    KPIs.Efficiency = W / (mfr_fuel * LHV);

    % Brake-Specific Fuel Consumption (BSFC)
    KPIs.BSFC = ((mfr_fuel*1000) * 3600) / (P / 1000); % g/kWh

    % Brake-Specific CO2 Emissions
    KPIs.bsCO2 = ((mfr_CO2*1000) * 3600) / (P / 1000); % g/kWh

    % Brake-Specific NOx Emissions
    KPIs.bsNOx = ((mfr_NOx*1000) * 3600) / (P / 1000); % g/kWh
end
