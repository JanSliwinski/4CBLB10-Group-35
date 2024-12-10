function [T_exh, Q_combustion_percycle, avg_m_fuelpercycle] = Texhaust(CA, C_p, mfr_fuel, mfr_air, RPM, W, LHV, T_int, Ncycles)

    % Inputs:
    % CA - Crank angle array
    % cp - 
    % mfr_fuel - Mass flow rate of fuel [g/s]
    % mfr_air - Mass flow rate of air [g/s]
    % RPM - Engine revolutions per minute
    % W - Work done [J]
    % LHV - Lower heating value of fuel [J/g]
    % T_int - Intake temperature [K]
    % Ncycles - Number of cycles
    
    % Outputs:
    % T_exh - Exhaust temperature array [K]
    
    
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
    
    %checking whether data is realistic
    avg_m_fuelpercycle = mean(m_fuel_percycle(2:end));
    avg_m_exh_percycle = mean(m_exh_percycle(2:end));
    avg_delt_T = mean(delta_T_percycle(2:end));
    avg_T_exh = mean(T_exh(2:end));
    avg_Q = mean(Q_combustion_percycle(2:end));
    disp(['Average fuel mass per cycle: ', num2str(avg_m_fuelpercycle)]);
   
    %disp(['Average exhaust mass per cycle: ', num2str(avg_m_exh_percycle)]);
    disp(['Average energy of combustion: ', num2str(avg_Q)]);
    disp(['Average delta T: ', num2str(avg_delt_T)]);
    disp(['Average temperature at exhaust: ', num2str(avg_T_exh)]);

    % Jings code for efficiency
    efficiency = (W / (mean(mfr_fuel(:)) * LHV))*100;
    disp(['Jing-s efficiency: ', num2str(efficiency), ' %']);

    % Plot the change of T_exhaust over all cycles
    figure;
    cycles = 1:Ncycles;
    plot(cycles, T_exh, 'LineWidth', 1);
    xlabel('Cycles');
    ylabel('Exhaust temperature (K)');
    xlim([1, 100]);
    title('Exhaust temperature over the 100 cycles run');
    grid on;
end
