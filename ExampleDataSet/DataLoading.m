% = ===========================================
% MATLAB Script: Data Integration and Loading
% ============================================
% This script loads experimental data TXT files and integrates additional
% measurements from a separate CSV file. It organizes the data based
% on Load (L), Intensity (I), and Composition (C) parameters for
% efficient analysis. Additionally, it applies a Savitzky-Golay filter
% to the Pressure data and stores the filtered results.
%
% Author: Ryan Sindic
% Date: December 20, 2024
% ============================================

%% Initialization
clear; clc; close all;

%% Define Paths
% Specify the folder containing the renamed experiment data TXT files
dataFolder = 'AdjustedData'; % <-- Replace with your actual folder path

% Specify the path to the additional data CSV file
additionalCSVPath = 'CompiledEmissions.csv'; % <-- Replace with your actual CSV file path

% Verify that the data folder exists
if ~isfolder(dataFolder)
    error('Data folder does not exist: %s', dataFolder);
end

% Verify that the additional CSV file exists
if ~isfile(additionalCSVPath)
    error('Additional data CSV file not found: %s', additionalCSVPath);
end

%% Define Savitzky-Golay Filter Parameters
% You can adjust these parameters as needed
k = 2;    % Polynomial order
n = 15;    % Window size (must be odd)
s = 0;    % Derivative order (0 for smoothing)

%% Load Additional Data
% Read the additional data CSV into a table
try
    additionalDataTable = readtable(additionalCSVPath, 'Delimiter', ';', 'ReadVariableNames', true);
    fprintf('Successfully loaded additional data CSV with %d rows.\n', height(additionalDataTable));
catch ME
    error('Failed to read additional data CSV: %s', ME.message);
end

% Validate that the required columns exist
requiredColumns = {'L', 'I', 'C', 'CO', 'HC', 'NOx', 'CO2', 'O2', 'Lambda'};
missingColumns = setdiff(requiredColumns, additionalDataTable.Properties.VariableNames);
if ~isempty(missingColumns)
    error('The additional CSV file is missing the following required columns: %s', strjoin(missingColumns, ', '));
end

% Handle missing Lambda values represented by '-'
if iscell(additionalDataTable.Lambda)
    additionalDataTable.Lambda(strcmp(additionalDataTable.Lambda, '-')) = {NaN};
    additionalDataTable.Lambda = cellfun(@str2double, additionalDataTable.Lambda);
elseif ischar(additionalDataTable.Lambda) || isstring(additionalDataTable.Lambda)
    additionalDataTable.Lambda = strrep(string(additionalDataTable.Lambda), '-', 'NaN');
    additionalDataTable.Lambda = str2double(additionalDataTable.Lambda);
end

% Create a unique key for each row in the additional data
additionalDataTable.UniqueKey = strcat('L', string(additionalDataTable.L), ...
                                      'I', string(additionalDataTable.I), ...
                                      'C', string(additionalDataTable.C));

% Initialize containers.Map for fast lookup
try
    additionalDataMap = containers.Map(additionalDataTable.UniqueKey, ...
                                       1:height(additionalDataTable));
    fprintf('Successfully created containers.Map for additional data.\n');
catch ME
    error('Failed to create containers.Map: %s', ME.message);
end

%% List Experiment Data Files
% Define the file pattern to match the experiment data TXT files
filePattern = fullfile(dataFolder, '*.txt'); % Adjust the extension if different

% Get a list of all matching files
files = dir(filePattern);

% Check if any files were found
if isempty(files)
    error('No experiment data TXT files found in folder: %s', dataFolder);
else
    fprintf('Found %d experiment TXT files in folder.\n', length(files));
end

%% Initialize Data Storage
% Initialize an empty cell array to store grouped data
% Columns:
% 1 - UniqueID (e.g., 'L2I25C30')
% 2 - Metadata [L, I, C]
% 3 - ExperimentData (cell array of matrices containing Crank Angle, Pressure, Current, Mass Flow)
% 4 - AdditionalData (struct)
dataArray = {}; 

%% Load Experiment Data with Parallel Processing
% Initialize a parallel pool if not already open
pool = gcp('nocreate'); % If no pool, do not create new one
if isempty(pool)
    parpool; % Opens the default parallel pool
    fprintf('Opened a new parallel pool.\n');
else
    fprintf('Using existing parallel pool.\n');
end

% Preallocate cell array for experiment data
experimentDataCell = cell(length(files), 1);

% Parallel loop to load data
parfor i = 1:length(files)
    fileName = fullfile(dataFolder, files(i).name);
    try
        % Read the experiment data
        % 'Delimiter',' ' for space-delimited, 'MultipleDelimsAsOne',true to handle multiple spaces
        tempData = readmatrix(fileName, ...
                              'Delimiter', ' ', ...
                              'MultipleDelimsAsOne', true);
        
        % Assign to cell array
        experimentDataCell{i} = tempData;
    catch ME
        warning('Failed to load file %s: %s', files(i).name, ME.message);
        experimentDataCell{i} = []; % Assign empty if failed to load
    end
end

%% Group Experiment Data
for i = 1:length(files)
    fileData = experimentDataCell{i};
    if isempty(fileData)
        fprintf('Skipping file %s: Failed to load data.\n', files(i).name);
        continue; % Skip files that failed to load
    end
    
    shortName = files(i).name;
    
    % Parse the File Name to Extract L, I, and C Values
    % Expected format: Ex{experiment}L{load}I{intensity}C{composition}.txt
    L_token = regexp(shortName, 'L(\d+(\.\d+)?)', 'tokens', 'once');
    I_token = regexp(shortName, 'I(\d+(\.\d+)?)', 'tokens', 'once');
    C_token = regexp(shortName, 'C(\d+(\.\d+)?)', 'tokens', 'once');
    
    % Ensure all tokens were found
    if isempty(L_token) || isempty(I_token) || isempty(C_token)
        fprintf('Skipping file %s: Filename does not match the expected pattern.\n', shortName);
        continue;
    end
    
    % Convert extracted values from cell to numeric
    L_value = round(str2double(L_token{1}));
    I_value = round(str2double(I_token{1}));
    C_value = round(str2double(C_token{1}));
    
    % Validate number of columns in fileData
    expectedCols = 4; % Crank Angle, Pressure, Current, Mass Flow
    if size(fileData, 2) ~= expectedCols
        fprintf('Skipping file %s: Expected %d columns, found %d columns.\n', shortName, expectedCols, size(fileData, 2));
        continue;
    end
    
    % Create a unique identifier based on L, I, and C
    uniqueID = sprintf('L%dI%dC%d', L_value, I_value, C_value);
    
    % Check if this uniqueID already exists in dataArray
    if isempty(dataArray)
        disp('before isempty')
        % dataArray is empty, add the first entry directly
        newRow = size(dataArray, 1) + 1;
        dataArray{newRow, 1} = uniqueID; % Store the unique identifier
        dataArray{newRow, 2} = [L_value, I_value, C_value]; % Store metadata
        dataArray{newRow, 3} = {fileData}; % Initialize column 3 as a cell array with the first entry
        dataArray{newRow, 4} = {}; % Initialize column 4 as empty, to be filled later
        fprintf('Creating new group %s (first entry).\n', uniqueID);
    else
         disp('after isempty')
        % dataArray is not empty, proceed to check for existing entries
        existingRow = find(strcmp({dataArray{:, 1}}, uniqueID));
        
        if ~isempty(existingRow)
            % If the group exists, append the data to the existing entry in column 3
            if iscell(dataArray{existingRow, 3})
                dataArray{existingRow, 3}{end+1, 1} = fileData;
                fprintf('Appending data to existing group %s.\n', uniqueID);
            else
                % Convert existing data to a cell array and add the new data
                dataArray{existingRow, 3} = {dataArray{existingRow, 3}, fileData};
                fprintf('Converting to cell array and appending data for group %s.\n', uniqueID);
            end
        else
            % If the group doesn't exist, create a new row
            newRow = size(dataArray, 1) + 1;
            dataArray{newRow, 1} = uniqueID; % Store the unique identifier
            dataArray{newRow, 2} = [L_value, I_value, C_value]; % Store metadata
            dataArray{newRow, 3} = {fileData}; % Initialize column 3 as a cell array with the first entry
            dataArray{newRow, 4} = {}; % Initialize column 4 as empty, to be filled later
            fprintf('Creating new group %s.\n', uniqueID);
        end
    end
end

%% Post-Grouping Verification
fprintf('After grouping, dataArray has %d rows and %d columns.\n', size(dataArray, 1), size(dataArray, 2));

% Optionally, display the first few entries
if ~isempty(dataArray)
    disp('First few entries in dataArray:');
    disp(dataArray(1:min(5, end), :));
else
    disp('dataArray is empty.');
end

%% Convert to MATLAB Table with Initial Data
% Define variable names for the table
variableNames = {'UniqueID', 'Metadata', 'ExperimentData', 'AdditionalData'};

% Check the size of dataArray
[numRows, numCols] = size(dataArray);
fprintf('Preparing to convert dataArray to table with %d rows and %d columns.\n', numRows, numCols);
fprintf('Number of VariableNames: %d\n', length(variableNames));

% Ensure dataArray has exactly 4 columns
if numCols ~= 4
    error('dataArray has %d columns, but 4 VariableNames were provided.', numCols);
end

% Convert to table
try
    T = cell2table(dataArray, 'VariableNames', variableNames);
    fprintf('Successfully converted dataArray to table T.\n');
catch ME
    error('Failed to convert cell array to table: %s', ME.message);
end

%% Assign Additional Data to Groups
for rowIdx = 1:size(dataArray, 1)
    currentL = dataArray{rowIdx, 2}(1);
    currentI = dataArray{rowIdx, 2}(2);
    currentC = dataArray{rowIdx, 2}(3);
    
    uniqueKey = sprintf('L%dI%dC%d', currentL, currentI, currentC);
    
    if isKey(additionalDataMap, uniqueKey)
        matchIdx = additionalDataMap(uniqueKey);
        
        % Handle multiple matches if necessary
        if length(matchIdx) > 1
            warning('Multiple matches found for group %s. Using the first match.', uniqueKey);
            matchIdx = matchIdx(1);
        end
        
        % Extract additional data columns
        CO = additionalDataTable.CO(matchIdx);
        HC = additionalDataTable.HC(matchIdx);
        NOx = additionalDataTable.NOx(matchIdx);
        CO2 = additionalDataTable.CO2(matchIdx);
        O2 = additionalDataTable.O2(matchIdx);
        Lambda = additionalDataTable.Lambda(matchIdx);
        
        % Store the additional data as a struct
        additionalData = struct('CO', CO, 'HC', HC, 'NOx', NOx, ...
                                'CO2', CO2, 'O2', O2, 'Lambda', Lambda);
        
        % Assign to column 4
        T.AdditionalData{rowIdx, 1} = additionalData;
    else
        % No matching additional data found
        warning('No additional data found for group %s (L=%.1f, I=%.1f, C=%.1f).', ...
                uniqueKey, currentL, currentI, currentC);
        T.AdditionalData{rowIdx, 1} = NaN; % Assign NaN or another placeholder
    end
end

%% Add FilteredPressure Column to Table
% Initialize the FilteredPressure column as a cell array
T.FilteredPressure = cell(height(T), 1);
fprintf('Added FilteredPressure column to the table.\n');

%% Apply Savitzky-Golay Filter to Pressure Data and Populate FilteredPressure Column
fprintf('Applying Savitzky-Golay filter to Pressure data...\n');
for rowIdx = 1:height(T)
    % Retrieve all ExperimentData matrices for the current group
    experiments = T.ExperimentData{rowIdx};
    
    % Initialize a cell array to store filtered pressure data for this group
    filteredPressureGroup = cell(size(experiments));
    
    for expIdx = 1:length(experiments)
        fileData = experiments{expIdx};
        
        % Extract Pressure data (assuming Pressure is the second column)
        pressureData = fileData(:, 2);
        
        % Apply the Savitzky-Golay filter
        try
            yyFilt = SGFilter(pressureData, k, n, s); % Using existing SGFilter function
            filteredPressureGroup{expIdx} = yyFilt;
        catch ME
            warning('Failed to filter Pressure data for group %s, experiment %d: %s', ...
                    T.UniqueID{rowIdx}, expIdx, ME.message);
            filteredPressureGroup{expIdx} = NaN; % Assign NaN if filtering fails
        end
    end
    
    % Assign the filtered Pressure data to the FilteredPressure column
    T.FilteredPressure{rowIdx, 1} = filteredPressureGroup;
    
    if mod(rowIdx, 10) == 0 || rowIdx == height(T)
        fprintf('Processed %d/%d groups.\n', rowIdx, height(T));
    end
end
fprintf('Completed applying Savitzky-Golay filter to all Pressure data.\n');

%% Define Crank Angle Resolution and Cycle Parameters
resolution = 0.2;  % Degrees crank angle resolution
NdatapointsPerCycle = 720 / resolution; % Number of data points per cycle

fprintf('Defined crank angle resolution: %.1f degrees.\n', resolution);
fprintf('Number of data points per cycle: %d.\n', NdatapointsPerCycle);

%% Add AverageCycleData Column to Table 
% Initialize the AverageCycleData column as a cell array
T.AverageCycleData = cell(height(T), 1);
fprintf('Added AverageCycleData column to the table.\n');

%% Compute Average Cycle Data for Each Group 
fprintf('Computing average cycle data for each group...\n');
for rowIdx = 1:height(T)
    % Retrieve all ExperimentData matrices and their corresponding filtered pressures
    experiments = T.ExperimentData{rowIdx};
    filteredPressures = T.FilteredPressure{rowIdx};
    
    % Initialize accumulators for sum
    sumPressure = zeros(NdatapointsPerCycle, 1);
    sumMassFlow = zeros(NdatapointsPerCycle, 1);
    sumCurrent = zeros(NdatapointsPerCycle, 1);
    cycleCount = 0;
    
    for expIdx = 1:length(experiments)
        fileData = experiments{expIdx};
        filteredPressure = filteredPressures{expIdx};
        
        % Determine the number of complete cycles in the data
        Nrows = size(fileData, 1);
        Ncycles = floor(Nrows / NdatapointsPerCycle);
        
        for cycle = 1:Ncycles
            startIdx = (cycle-1)*NdatapointsPerCycle + 1;
            endIdx = cycle * NdatapointsPerCycle;
            
            % Extract data for the current cycle
            try
                pressureCycle = filteredPressure(startIdx:endIdx);
                massFlowCycle = fileData(startIdx:endIdx, 4);
                currentCycle = fileData(startIdx:endIdx, 3);
                
                % Accumulate the data
                sumPressure = sumPressure + pressureCycle;
                sumMassFlow = sumMassFlow + massFlowCycle;
                sumCurrent = sumCurrent + currentCycle;
                
                cycleCount = cycleCount + 1;
            catch ME
                warning('Failed to process cycle %d in group %s, experiment %d: %s', ...
                        cycle, T.UniqueID{rowIdx}, expIdx, ME.message);
                continue; % Skip this cycle if there's an error
            end
        end
    end
    
    if cycleCount > 0
        % Compute the average values
        avgPressure = sumPressure / cycleCount;
        avgMassFlow = sumMassFlow / cycleCount;
        avgCurrent = sumCurrent / cycleCount;
        
        % Store the averaged data in a struct
        averageData = struct('AvgPressure', avgPressure, ...
                             'AvgMassFlow', avgMassFlow, ...
                             'AvgCurrent', avgCurrent);
    else
        % If no cycles were processed, assign NaN
        averageData = NaN;
        warning('No complete cycles found for group %s. AverageCycleData set to NaN.', T.UniqueID{rowIdx});
    end
    
    % Assign the averaged data to the new column
    T.AverageCycleData{rowIdx, 1} = averageData;
    
    if mod(rowIdx, 10) == 0 || rowIdx == height(T)
        fprintf('Processed %d/%d groups.\n', rowIdx, height(T));
    end
end
fprintf('Completed computing average cycle data for all groups.\n');


%% Display the Table
disp('Integrated Data Table:');
disp(T(1:min(5, height(T)), :)); % Display first 5 rows


%% Cleanup
% Optionally, shut down the parallel pool if no longer needed
% delete(gcp('nocreate'));

% ============================================
% End of Script
% ============================================
