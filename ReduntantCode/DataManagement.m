% Specify location of the experimental file(s) and open them
folderPath = './Data/session1_Raw';
outputfilePath = './Data/Processed_session1';

n_cycles = 100;
cycle_length = 3600;

% If averaging data also, use this code:
% sumData = 0;

% Create output directory if it doesn't exist
if ~isfolder(outputfilePath)
    mkdir(outputfilePath);
end

% Get a list of all .txt files in the folder
fileList = dir(fullfile(folderPath, '*.txt'));
numFiles = length(fileList);

if isempty(fileList)
    error('No .txt files found in the folder: %s', folderPath);
end

% Loop through each file and accumulate data
for magic1 = 1:numFiles
    try
        filePath = fullfile(folderPath, fileList(magic1).name);
        disp(['Processing file: ', filePath]);

        % Load data
        data = readmatrix(filePath);

        % Ensure sufficient columns exist
        if size(data, 2) < 4
            error('File %s does not have enough columns.', filePath);
        end

        % Apply Savitzky-Golay filter
        for magic2 = 2:4
            column = data(:,magic2);
            filtered = sgolayfilt(column, 2, 15); 
            data(:,magic2) = filtered;
        end

        % Adjust per cycle
        for magic3 = 0:99 
            startIDX = 3600*magic3+1; 
            endIDX = 3600*(magic3+1);
            midIDX = round((startIDX + endIDX) / 2);

            if data(midIDX,2) < 10
                data(startIDX:endIDX,:) = NaN;
            else
                bias = min(data(startIDX:endIDX,2));
                data(startIDX:endIDX,2) = data(startIDX:endIDX,2) - bias + 1;
            end
        end

        % Remove rows set to NaN
        data(any(isnan(data), 2), :) = [];

        % Using a vectorized approach to average cycles
        data_reshaped = reshape(data, cycle_length, n_cycles, 4);
        data_avg = mean(data_reshaped, 2);
        data = squeeze(data_avg);


        % Save filtered data
        [~, name, ext] = fileparts(fileList(magic1).name);
        outputFileName = fullfile(outputfilePath, [name, '_filt_avg', ext]);
        disp(['Saving to: ', outputFileName]);
        writematrix(data, outputFileName);
    catch ME
        warning('Error processing file %s: %s', fileList(magic1).name, ME.message);
    end
    % sumData = sumData + data;
end


% averageData = sumData / numFiles;

% Save averaged data
% outputFileName = fullfile(outputfilePath, 'averaged_filtered_data_2_IMEP.txt');
% disp(['Saving averaged data to: ', outputFileName]);
% writematrix(averageData, outputFileName);