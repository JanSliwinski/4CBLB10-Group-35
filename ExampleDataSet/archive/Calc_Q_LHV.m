function [T_exh, Q_combustion_LHV, m_combusted] = Texhaust(C_p, mfr_fuel, mfr_air, RPM, W, LHV, T_int)
    % Inputs:
    % CA - Crank angle array
    % cp - specific heat of exhaust gasses as calculated in cp_cv_gamma.m
    % mfr_fuel - Mass flow rate of fuel [g/s]
    % mfr_air - Mass flow rate of air [g/s]
    % RPM - Engine revolutions per minute
    % W - Work done [J]
    % LHV - Lower heating value of fuel [J/g]
    % T_int - Intake temperature [K]
    % Ncycles - Number of cycles
    
    % Outputs:
    % T_exh - Exhaust temperature array [K]
    % Q_combustion_LHV - Energy of combustion trhough the LHV way
    % m_combusted - mass of fuel that takes part in the combustion
    
    % Calculate time step per cycle
    t_cycle = (60/RPM)*2; % time per cycle for a four stroke engine [s]
    
    % Calculate per-cycle fuel and exhaust masses
    m_combusted = mfr_fuel * t_cycle;   %mass of fuel combusted in a cycle [g]
    mfr_exh = mfr_air + mfr_fuel;   %total exhaust mass flow rate [g/s]
    m_exh = mfr_exh * t_cycle;      %mass of exhaust in a cycle [g]

    % Calculate combustion energy and heat for warming exhaust
    Q_combustion_LHV = m_combusted * LHV; % Combustion energy per cycle [J]
    Q_warmingtheair = Q_combustion_LHV - W; % Energy for warming exhaust gases [J]
    
    % Compute temperature change and exhaust temperature
    delta_T_percycle = Q_warmingtheair / (C_p * m_exh); % Temp change per cycle [K]
    T_exh = T_int + delta_T_percycle; % Exhaust temperature per cycle [K]
    
    %Display data for ease of checking:
    disp(['Fuel mass combusted per cycle: ', num2str(m_combusted), 'g']);
    disp(['Duration of one cycle: ', num2str(t_cycle), 's']);
    disp(['Energy of combustion-LHV: ', num2str(Q_combustion_LHV), 'J']);
    disp(['Delta T - LHV: ', num2str(delta_T_percycle), 'K']);
    disp(['Temperature at exhaust - LHV: ', num2str(T_exh), 'K']);
end
