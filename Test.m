clear all;
close all;

%xData = getData('xData.xlsx');


%TetheredParticleAnalysis_Original(xData, 200, 100 , 0.05);
%figure;
%TetheredParticleAnalysis(xData, 200, 100 , 0.04);
%figure;
%TetheredParticleAnalysis_Version2(xData, 200, 200, 1e-12 ,[0, 1], 0.01);
%figure;
%TetheredParticleAnalysis_Version2(xData, 200, 200, 0.5e-12 ,[0.3, 0.65], 0.01);
%figure;
%f <= 25N ~ deltal = 0.001
laserForce = 15e-12; %Newton
appliedRange = [0.2, 0.86];
%TetheredParticleAnalysis_Version2(xData, 500, 150, laserForce ,appliedRange, 0.0001);
%TetheredParticleAnalysis_WithDependentTimestep(200, 100, laserForce,pi/4, appliedRange, 0.0000001, 50000*2);
TetheredParticleAnalysis_WithDependentTimestep_3D(200, 100, laserForce,0, 0, appliedRange, 0.001, 35000);
%The problem that makes extension exceeds Lo is that we chose too large
%time step
%Becuase of that, it's not continous enough to simulate the extendsion
%between step
%Cause when force is much larger, it accelerate the bead faster, we need
%smaller time step, which means the system must capture faster
%There must be some correlation between applied force and time step