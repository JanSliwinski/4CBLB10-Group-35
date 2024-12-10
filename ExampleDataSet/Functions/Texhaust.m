function [T_exh, avg_metrics, efficiency] = Texhaust(CA, mfr_fuel, mfr_air, RPM, W, LHV, T_int, Ncycles)
    % Texhaust calculates the exhaust temperature and related metrics
    % Inputs:
    % CA - Crank angle array
    % mfr_fuel - Mass flow rate of fuel [g/s]
    % mfr_air - Mass flow rate of air [g/s]
    % RPM - Engine revolutions per minute
    % W - Work done [J]
    % LHV - Lower heating value of fuel [J/g]
    % T_int - Intake temperature [K]
    % Ncycles - Number of cycles
    
    % Outputs:
    % T_exh - Exhaust temperature array [K]
    % avg_metrics - Struct containing average metrics
    % efficiency - Efficiency of the engine process [%]
    
    % Constants
    C_p = 1.005; % Specific heat capacity of air [J/g*K]
    
    % Calculate time step per crank angle
    DeltaCA = mean(diff(CA(:,1))); % Step size between crank angles (degrees)
    Delta_t_perCA = (60 * DeltaCA) / (RPM * 360); % Time per crank angle step (s)
    
    % Calculate per-cycle fuel and exhaust masses
    m_fuel_percycle = sum(mfr_fuel, 1) * Delta_t_perCA; % Fuel mass per cycle [g]
    mfr_exh = mfr_fuel + mfr_air; % Exhaust mass flow rate [g/s]
    m_exh_percycle = sum(mfr_exh, 1) * Delta_t_perCA; % Exhaust mass per cycle [g]
    
    % Calculate combustion energy and heat for warming exhaust
    Q_combustion_percycle = m_fuel_percycle * LHV; % Combustion energy per cycle [J]
    Q_warmingtheair = Q_combustion_percycle - W; % Energy for warming exhaust gases [J]
    
    % Compute temperature change and exhaust temperature
    delta_T_percycle = Q_warmingtheair ./ (C_p * m_exh_percycle); % Temp change [K]
    T_exh = T_int + delta_T_percycle; % Exhaust temperature per cycle [K]
    
    % Average metrics calculation
    avg_metrics.avg_m_fuel_percycle = mean(m_fuel_percycle(2:end));
    avg_metrics.avg_m_exh_percycle = mean(m_exh_percycle(2:end));
    avg_metrics.avg_delta_T = mean(delta_T_percycle(2:end));
    avg_metrics.avg_T_exh = mean(T_exh(2:end));
    avg_metrics.avg_Q_combustion = mean(Q_combustion_percycle(2:end));
    
    % Efficiency calculation
    efficiency = (W / (mean(mfr_fuel(:)) * LHV)) * 100;
    
    % Display metrics
    disp(['Average fuel mass per cycle: ', num2str(avg_metrics.avg_m_fuel_percycle)]);
    disp(['Average exhaust mass per cycle: ', num2str(avg_metrics.avg_m_exh_percycle)]);
    disp(['Average combustion energy: ', num2str(avg_metrics.avg_Q_combustion)]);
    disp(['Average delta T: ', num2str(avg_metrics.avg_delta_T)]);
    disp(['Average exhaust temperature: ', num2str(avg_metrics.avg_T_exh)]);
    disp(['Efficiency: ', num2str(efficiency), ' %']);
end
