clear all, close all;
load("T.mat");
gamma = CalculateGamma(T,3);
figure;plot(0:0.2:720-0.2,gamma);