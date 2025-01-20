function [SpS, El] = myload(filename, species)
%MYLOAD Loads NASA species and element data from a MAT-file for specified species.
%
% Inputs:
%   filename : Path to the MAT-file containing variables 'Sp' and 'El'
%   species  : Cell array of species names to select
%
% Outputs:
%   SpS      : Array of species structures corresponding to specified species
%   El       : Element structure array with masses converted to kilograms
%
% This function loads species (Sp) and element (El) data from the given file,
% selects specific species based on input, and adjusts element masses for consistency.

    %% Load Data from File
    % Load the MAT-file which contains variables 'Sp' and 'El'
    load(filename);

    %% Select Specified Species
    % Find indices of species in 'Sp' that match the provided list of species names
    isp = myfind({Sp.Name}, species);

    % Extract species structures corresponding to the found indices
    SpS = Sp(isp);

    %% Adjust Element Masses
    % Convert element masses from grams to kilograms for each element in 'El'
    for i = 1:length(El)
        El(i).Mass = El(i).Mass * 1e-3;
    end
end
