warning off
clc, clear all, close all;

%% Loading Data
DataArray = cell(0,3);
folderPath = 'Data'; 
filePattern = fullfile(folderPath, '*.txt'); 
files = dir(filePattern);
% Loop through each file in the folder
for i = 1:length(files)
    % Get the full path of the current file
    fileName = fullfile(folderPath, files(i).name);
    
    dataFileName = fullfile('Data' , 'processed_Data_experiment1_load3.5.txt');
    dataIn = table2array(readtable(dataFileName));

    %Error handling:
    [Nrows, Ncols] = size(dataIn);
    if Ncols  ~= 4
        warning('Data loaded does not have the expected amount of columns (4)');
    end

    resolution = 0.2;  % Degrees crank angle resolution
    NdatapointsPerCycle = 720 / resolution;
    Ncycles = Nrows / NdatapointsPerCycle;

    if mod(Nrows, NdatapointsPerCycle) ~= 0
        error('Number of data points is not an integer multiple of data points per cycle.');
    end

    % Reshape data into cycles
    Ca = reshape(dataIn(:, 1), [], Ncycles);      % Crank angle in degrees
    p = reshape(dataIn(:, 2), [], Ncycles) * 1e5;  % Pressure in Pa
    S_current = reshape(dataIn(:, 3), [], Ncycles);  % Sensor current 
    mfr_fuel = reshape(dataIn(:, 4), [], Ncycles);  % Fuel mass flow
    fileData = [Ca , p , S_current, mfr_fuel];
    % Extract the file name without the path
    shortName = files(i).name;
    
    % **Parse the File Name to Extract L, I, and C Values**
    % Assuming the format is Ex{num}L{num}I{num}C{num}, e.g., Ex1L50I21C30
    L_value = regexp(shortName, 'L(\d+)', 'tokens', 'once');
    I_value = regexp(shortName, 'I(\d+)', 'tokens', 'once');
    C_value = regexp(shortName, 'C(\d+)', 'tokens', 'once');
    
    % Convert extracted values from cell to numeric
    L_value = str2double(L_value{1});
    I_value = str2double(I_value{1});
    C_value = str2double(C_value{1});
    
    % Create a unique identifier based on L, I, and C
    uniqueID = sprintf('L%dI%dC%d', L_value, I_value, C_value);
    
    % **Check if this uniqueID already exists in dataArray**
    existingRow = find(strcmp({dataArray{:, 1}}, uniqueID));
    
    if ~isempty(existingRow)
        % If the group exists, append the data to the existing entry in column 3
        % Initialize as a cell array if it's the first addition
        if iscell(dataArray{existingRow, 3})
            dataArray{existingRow, 3}{end+1, 1} = fileData;
        else
            % Convert existing data to a cell array and add the new data
            dataArray{existingRow, 3} = {dataArray{existingRow, 3}, fileData};
        end
    else
        % If the group doesn't exist, create a new row
        newRow = size(dataArray, 1) + 1;
        dataArray{newRow, 1} = uniqueID; % Store the unique identifier
        dataArray{newRow, 2} = [L_value, I_value, C_value]; % Optional: Store metadata
        dataArray{newRow, 3} = {fileData}; % Initialize column 3 as a cell array with the first entry
        dataArray{newRow, 4} = {}; % Initialize column 4 as empty
    end
end


additionalCSVPath = 'path_to__additional_data.csv'; 
additionalDataTable = readtable(additionalCSVPath);

% Validate columns
requiredColumns = {'L', 'I', 'C', 'CO', 'HC', 'NOx', 'CO2', 'O2', 'Lambda'};
missingColumns = setdiff(requiredColumns, additionalDataTable.Properties.VariableNames);
if ~isempty(missingColumns)
    error('The additional CSV file is missing the following required columns: %s', strjoin(missingColumns, ', '));
end

for rowIdx = 1:size(dataArray, 1)
    currentL = dataArray{rowIdx, 2}(1);
    currentI = dataArray{rowIdx, 2}(2);
    currentC = dataArray{rowIdx, 2}(3);
    
    matchIdx = find(additionalDataTable.L == currentL & ...
                    additionalDataTable.I == currentI & ...
                    additionalDataTable.C == currentC);
    
    if ~isempty(matchIdx)
        % Handle multiple matches if necessary
        if length(matchIdx) > 1
            warning('Multiple matches found for group %s (L=%d, I=%d, C=%d). Using the first match.', ...
                    dataArray{rowIdx, 1}, currentL, currentI, currentC);
            matchIdx = matchIdx(1); % Use the first match
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
        dataArray{rowIdx, 4} = additionalData;
    else
        warning('No additional data found for group %s (L=%d, I=%d, C=%d).', ...
                dataArray{rowIdx, 1}, currentL, currentI, currentC);
        dataArray{rowIdx, 4} = NaN; % Or another placeholder
    end
end
