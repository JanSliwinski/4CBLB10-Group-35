%% Function to compute cp based on T
function cv = compute_cv(T, SpS, massComposition)
% COMPUTE_CV Calculates the specific heat capacity at constant volume for a mixture
%
% This function computes the specific heat capacity at constant volume (cv) 
% for a given mixture of species at a specified temperature using NASA 
% thermodynamic database.
%
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
%   cv              - Specific heat capacity at constant volume 
%                     Type: Scalar numeric value
%                     Units: J/(kg·K)
%                     Calculation: Weighted average of individual species cv
%
% Method:
%   1. Loads NASA thermal database containing species-specific thermodynamic data
%   2. Calculates individual species cv using CvNasa function
%   3. Converts individual species cv from kJ/(mol·K) to J/(kg·K)
%   4. Computes mixture cv as a mass-weighted sum of individual species cv
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
%   result_cv = compute_cv(T, SpS, massComp);

%% Load Nasa database
TdataBase = fullfile('Nasa', 'NasaThermalDatabase');
load(TdataBase);


    NSp = length(SpS);
    cvi = zeros(NSp, 1);
    for i = 1:NSp
        % CpNasa returns cp in kJ/(mol·K), convert to J/(kg·K)
        cvi(i) = (CvNasa(T, SpS(i)));  
    end
    % Compute mixture cp as weighted sum of species cp
    cv = massComposition * cvi;  % Result in J/(kg·K)
end