clear all;
close all;
%% Params: 

laserForce = 5e-12; %Newton
appliedRange = [0.33, 0.66];
alpha = 0; %theta
beta = 0; %phi

Lo = 40000;
Rb = 1500;

eta = 2*10^(-30);
temperature = 25; %Celcius DE
TetheredParticleAnalysis_WithDependentTimestep_3D(laserForce,alpha ,beta , appliedRange, 0.01, 8000, Lo, Rb, eta, temperature);
