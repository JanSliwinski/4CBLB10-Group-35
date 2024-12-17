function true_mass_flow = CalculateMassFlowFuel(mfr_fuel, S_current, Ca, RPM, threshold)
%% Function to calculate the true mass flow rate during the injection period
%
% Inputs:
%   - mfr_fuel: Matrix of fuel mass flow rates (rows = data points, columns = cycles)
%   - S_current: Matrix of sensor current (rows = data points, columns = cycles)
%   - Ca: Matrix of crank angles (rows = data points, columns = cycles)
%   - RPM: Engine speed in revolutions per minute
%   - threshold: Sensor current threshold to detect injection start and end
%
% Output:
%   - true_mass_flow: True mass flow rate during the injection period (g/s)

    % Number of cycles
    Ncycles = size(mfr_fuel, 2);

    % Initialize an array to store the true mass flow rates for each cycle
    true_mass_flow_per_cycle = zeros(Ncycles, 1);

    % Calculate the time per crank angle step (in seconds)
    time_per_cycle = 2 * (60 / RPM);  % Time for one full engine cycle (720°) in seconds
    crank_angle_step = Ca(2, 1) - Ca(1, 1);  % Crank angle step in degrees
    time_per_step = time_per_cycle / (720 / crank_angle_step);  % Time per step in seconds

    for i = 1:Ncycles
        % Extract sensor current and fuel mass flow for the current cycle
        S_current_cycle = S_current(:, i);
        mfr_fuel_cycle = mfr_fuel(:, i);

        % Detect the injection start and end indices where the sensor current exceeds the threshold
        injection_indices = find(S_current_cycle > threshold);

        if isempty(injection_indices)
            warning('No injection detected in cycle %d', i);
            continue;
        end

        injection_start_idx = injection_indices(1);
        injection_end_idx = injection_indices(end);
        

        % Calculate the injection duration in seconds
        injection_duration = (injection_end_idx - injection_start_idx + 1) * time_per_step;

        % Calculate the total mass injected per cycle (assuming mfr_fuel is in g/s)
        total_mass_injected = sum(mfr_fuel_cycle(injection_start_idx:injection_end_idx)) * time_per_step; % Total mass injected in g

        % Calculate the true mass flow rate during injection
        true_mass_flow_per_cycle(i) = total_mass_injected / injection_duration;
    end

    % Calculate the overall average true mass flow rate across all cycles
    true_mass_flow = mean(true_mass_flow_per_cycle, 'omitnan');

    % Display the result
    fprintf('True mass flow rate during injection: %.6f g/s\n', true_mass_flow);
end
