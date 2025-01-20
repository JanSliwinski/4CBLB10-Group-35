clear all, close all; addpath("Functions\")
load("T.mat"); resolution = 0.2;

        %% Load Nasa database
        addpath("Nasa\", "Functions\");
        TdataBase = fullfile('Nasa', 'NasaThermalDatabase'); 
        load(TdataBase);
        
        %% Find species
        iSp = myfind({Sp.Name}, {'N2', 'O2', 'CO2','H2O','Diesel'}); % Find indexes of these species
        SpS = Sp(iSp); % Create subset of species based on found indexes
        
gamma = CalculateGamma(T,3,SpS);

%%

aROHR = get_aROHR(T,3,gamma);

% Crop aROHR
aROHR = aROHR(345/resolution:405/resolution);
aHR = get_aHR(aROHR);

figure;
subplot(1,2,1)
plot(linspace(345,405,301),aROHR)
subplot(1,2,2)
plot(linspace(345,405,301),aHR)
