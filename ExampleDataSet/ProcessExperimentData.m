% function to process the experimental data by averaging every single entry
% and outputting a .txt file containing 4 columns that can be imediately
% used in the main code

clc, clear all, close all

% Folder containing the .txt files the data from the first experiment
folderPath = './ExampleDataSet/Data/session1_Raw'; %defines folder path relative to the current directory

% Get a list of all .txt files in the folder
fileList = dir(fullfile(folderPath, '*.txt'));
numFiles = length(fileList); %number of files

%Read the first file to see how many rows our output file must have
% (assume that all .txt raw files have the same amount of rows and columns)
FirstFilePath = fullfile(folderPath, fileList(1).name);
sampleData = load(FirstFilePath);
[numRows, numCols] = size(sampleData);

%initialise a dataset with this size that will store as it's entries the
%sum of all corresponding (17) entries of all data files
sumData = zeros(numRows, numCols);

% Loop through each file and accumulate data
for i = 1:numFiles
    filePath = fullfile(folderPath, fileList(i).name); % Construct full file path to each file

    data = load(filePath);   % Read the data from the file
    
    % Check if the data has the same dimensions as the sample file
    if size(data, 1) ~= numRows || size(data, 2) ~= numCols
        error('File %s dimensions do not match other files.', fileList(i).name);
    end
    
    % Accumulate data
    sumData = sumData + data;
end

%Average each entry knowing the number of exprtimental data files
averageData = sumData / numFiles;

%Write this avaraged data to a new .txt file so it can be easily
%implemented in the main code
writematrix(averageData, './ExampleDataSet/Data/processed_Data_experiment1.txt');

%NOTE: the columns contain the following data:
% 1 - Crank Angle (deg)
% 2 - Pressure (bar)
% 3 - Sensor current (mA)
% 4 - Fuel mass flow (g/s)

% Display success message
fprintf('Processed data saved to: %s\n', outputFilePath);


