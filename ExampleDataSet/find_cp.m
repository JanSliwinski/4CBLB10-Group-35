clc, clear all, close all; 

%% Add necessary paths
relativepath_to_generalfolder = 'Nasa'; % Adjust if necessary
addpath(relativepath_to_generalfolder);

%% Load Nasa database
TdataBase = fullfile('Nasa', 'NasaThermalDatabase');
load(TdataBase);

%% Find species
iSp = myfind({Sp.Name}, {'O2', 'CO2', 'N2'});  % Find indexes of these species
SpS = Sp(iSp);  % Subset of species
NSp = length(SpS);

% Given variables
LHV = 50 * 1e6;         % Lower Heating Value in J/kg
m_dot = 0.00013;        % Mass flow rate (kg/s)
Q_dot = LHV * m_dot;    % Heat transfer rate (W)
T_initial = 295.15;     % Initial temperature (K)
tolerance = 1000;          % Acceptable error in W
iterations = 100000;
deltaT = 1;           % Initial guess for deltaT (K)
error = Inf;

% Mole fractions (assumed to be given)
O2_percent = 0.1442;
CO2_percent = 0.04667;
N2_percent = 1 - O2_percent - CO2_percent;
moleFractions = [O2_percent, CO2_percent, N2_percent];

% Convert mole fractions to mass fractions
Mi = [SpS.Mass] / 1000;  % Molar masses in kg/mol
massComposition = (moleFractions .* Mi) / sum(moleFractions .* Mi);

counter = 0;
while error > tolerance
    T = T_initial + deltaT;
    
    % Compute cp based on average T (in J/kg·K)
    T_avg = (T + T_initial) / 2;
    cp = compute_cp(T_avg, SpS, massComposition);
    
    % Calculate Q̇ based on current cp and deltaT
    Q_dot_calculated = m_dot * cp * deltaT;
    
    % Compute the error between calculated and given Q̇
    error = abs(Q_dot_calculated - Q_dot);
    
    % Increment deltaT
    deltaT = deltaT + 10;
    counter = counter + 1;
end

disp("Great Success, cp = ");
disp(cp);

% Find the deltaT with minimum error
% % % [min_error, idx] = min(error);
% % % if min_error <= tolerance
% % %     deltaT_solution = delta(idx);
% % %     T_solution = T_initial + deltaT_solution;
% % %     cp_solution = compute_cp((T_initial + T_solution)/2, SpS, massComposition);
% % %     fprintf('Found solution after %d iterations.\n', idx);
% % %     fprintf('deltaT = %.5f K\n', deltaT_solution);
% % %     fprintf('Final T = %.5f K\n', T_solution);
% % %     fprintf('cp = %.5f J/kg·K\n', cp_solution);
% % % else
% % %     fprintf('Solution not found within maximum iterations.\n');
% % % end


%% Function to compute cp based on T
function cp = compute_cp(T, SpS, massComposition)
    NSp = length(SpS);
    Mi = [SpS.Mass] / 1000;  % Molar masses in kg/mol
    cpi = zeros(NSp, 1);
    for i = 1:NSp
        % CpNasa returns cp in J/(mol·K), convert to J/(kg·K)
        cpi(i) = CpNasa(T, SpS(i))/Mi(i);
    end
    % Compute mixture cp as weighted sum of species cp
    cp = 0.001 .* massComposition * cpi;  % Result in J/kg·K
end
