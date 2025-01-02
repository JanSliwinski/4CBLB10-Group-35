clear all, close all; addpath("Functions\")
load("T.mat"); resolution = 0.2;
gamma = CalculateGamma(T,3);

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
