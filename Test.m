clear all;
close all;
%% Params: 

laserForce = 5e-12; %Newton
appliedRange = [0.05, 0.95];
alpha = 0; %theta
beta = 0; %phi

Lo = 1200;
Rb = 120;

eta = 2*10^(-30);
temperature = 25; %Celcius DE
dataPoints = 5000;
TetheredParticleAnalysis_WithDependentTimestep_3D(laserForce,alpha ,beta , appliedRange, 0.01, dataPoints, Lo, Rb, eta, temperature);
