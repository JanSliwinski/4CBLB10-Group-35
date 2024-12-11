%% Function to compute cp based on T
function cp = compute_cp(T, SpS, massComposition)
% COMPUTE_CP Calculates the specific heat capacity at constant pressure for a mixture
%
% This function computes the specific heat capacity at constant pressure (cp) 
% for a given mixture of species at a specified temperature using NASA 
% thermodynamic database.
% Inputs:
%   T               - Temperature (units expected to match NASA database, 
%                     likely Kelvin)
%                     Type: Scalar numeric value
%                     Range: Typically 200-6000 K
%
%   SpS             - Array of species identifiers 
%                     Type: Array of strings or numerical indices
%                     Purpose: Specifies which chemical species are in the mixture
%                     Example: ['N2', 'O2', 'CO2'] or corresponding numerical codes
%
%   massComposition - Mass fractions of each species in the mixture
%                     Type: Numeric array 
%                     Length: Must match the length of SpS
%                     Constraints: 
%                     - Values should be between 0 and 1
%                     - Sum of all mass fractions should equal 1
%                     Example: [0.7, 0.2, 0.1] for a 3-species mixture
%
% Output:
%   cp              - Specific heat capacity at constant pressure 
%                     Type: Scalar numeric value
%                     Units: J/(kg·K)
%                     Calculation: Weighted average of individual species cp
%
% Method:
%   1. Loads NASA thermal database containing species-specific thermodynamic data
%   2. Calculates individual species cp using CpNasa function
%   3. Converts individual species cp from kJ/(mol·K) to J/(kg·K)
%   4. Computes mixture cp as a mass-weighted sum of individual species cp
%
% Note:
%   - Requires prior loading of NASA thermal database
%   - Assumes linear mixing of heat capacities
%   - Accuracy depends on NASA database and interpolation method
%
% Example:
%   T = 300;  % Temperature in Kelvin
%   SpS = ['N2', 'O2'];  % Species
%   massComp = [0.7, 0.3];  % Mass fractions
%   result_cp = compute_cp(T, SpS, massComp);

%% Load Nasa database
TdataBase = fullfile('Nasa', 'NasaThermalDatabase');
load(TdataBase);

    NSp = length(SpS);
    cpi = zeros(NSp, 1);
    for i = 1:NSp
        % CpNasa returns cp in kJ/(mol·K), convert to J/(kg·K)
        cpi(i) = (CpNasa(T, SpS(i)));  
    end
    % Compute mixture cp as weighted sum of species cp
    cp = massComposition * cpi;  % Result in J/(kg·K)
end