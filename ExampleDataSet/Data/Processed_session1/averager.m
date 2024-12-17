clear all

% Get a list of all .txt files in the folder
fileList = dir(fullfile('*.txt'));
numFiles = length(fileList);

% Loop through each file and accumulate data
for magic1 = 1:numFiles
    try
        filePath = fullfile(fileList(magic1).name);
        disp(['Processing file: ', filePath]);

        % Load data
        data = readmatrix(filePath);

        % Using a vectorized approach to average cycles
        n_cycles = 100;
        cycle_length = 3600;
        data_reshaped = reshape(data, cycle_length, n_cycles, 4);
        data_avg = mean(data_reshaped, 2);
        data = squeeze(data_avg);

        % Save filtered data
        [~, name, ext] = fileparts(fileList(magic1).name);
        outputFileName = fullfile([name, '_averaged', ext]);
        disp(['Saving to: ', outputFileName]);
        writematrix(data, outputFileName);
    catch ME
        warning('Error processing file %s: %s', fileList(magic1).name, ME.message);
    end
    % sumData = sumData + data;
end
