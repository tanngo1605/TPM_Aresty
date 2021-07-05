
function TetheredParticleAnalysis_WithDependentTimestep_3D (fLaser,alpha,beta,appliedRange,deltat,datapoints, Lo, Rb, eta, temperature)
%theta: angle with z axis
%phi: angle in x0y plane
%alpha: angle of force with z axis
%beta: angle of the force in xOy plane.

n         = datapoints;    %number of data points
%This line of code might not be needed
time = (1:n)*deltat;

%%%%Simulated data%%%%
[Xsim,Ysim,Zsim, FullExtension, AngleThetaSimDegree, AnglePhiSimDegree, DNAForce, timeline] = RandWalkSim(n, fLaser,alpha, beta, appliedRange, Lo, Rb, eta, temperature);  %generate simulated data


%%%%output%%%%
time = timeline;%time series
subplot(3,3,1);
plot(time,Xsim,'b')%Experimental data
title ('X coordinate data');
xlabel ('Time (sec)');
ylabel ('X Bead Position (nm)');

subplot(3,3,2);
plot(time, Ysim,'r')%Experimental data
title ('Y coordinate data');
xlabel ('Time (sec)');
ylabel ('Y Bead Position (nm)');

subplot(3,3,3);
plot(time, Zsim,'b')%Experimental data
title ('Z coordinate data');
xlabel ('Time (sec)');
ylabel ('Z Bead Position (nm)');


subplot(3,3,4);
plot(time, FullExtension,'r')%Experimental data
title ('Full extension');
xlabel ('Time (sec)');
ylabel ('Extension Bead Position (nm)');

subplot(3,3,5);
plot(time, AngleThetaSimDegree,'b')%Experimental data
title ('Theta degree');
xlabel ('Time (sec)');
ylabel ('Degree (o)');

subplot(3,3,6);
plot(time, AnglePhiSimDegree,'b')%Experimental data
title ('Phi degree');
xlabel ('Time (sec)');
ylabel ('Degree (o)');

subplot(3,3,7);
plot(time, DNAForce,'b')%Experimental data
title ('Stretching DNA force');
xlabel ('Time (sec)');
ylabel ('Magnitude (N)');


function [Xsim, Ysim, Zsim, FullExtension, AngleSimThetaDegree , AngleSimPhiDegree, DNAForce, timeline]=RandWalkSim(n,fLaser,alpha, beta,appliedRange, Lo, Rb, eta, temperature)
%This function simulates a 1D random walk in a harmonic potential
%%%%physical parameters%%%%
% Lo       = 3477*.34;            % Lo = tether length (nm) = 1182.18
Lp      = 72;                   % Lp = tether persistence length (nm)
% Rb      = 240;                  %bead radius (nm)
% Rb = 50;
%change to dynamic kbT = k * T
% kbT     = 4.1*10^(-21);         %thermal energy (J)
kbT = (1.38*10^-23) * (temperature + 273);
% eta     = 2.4*10^(-30);         %viscosity of H2O (J*s/nm^3)
D       = kbT/(6*pi*eta*Rb);    %Stokes diffusion constant (nm2/s)
sita0    = 6*pi*eta*Rb;          %drag coefficient, correct for close proximity for surface. FAXEN's law, 
deltat = getNextTimestep(sita0, Lp, Lo, 0, kbT);
%deltat = 0.001;
Xsim=zeros(1,n);
Ysim=zeros(1,n);
Zsim=zeros(1,n);
Zsim =  Zsim + Rb;

DNAForce = zeros(1,n);
FullExtension = zeros(1,n);
AnglePhiSim = zeros(1,n);
AngleThetaSim = zeros(1,n);
AngleThetaSim = AngleThetaSim + pi/2;
AngleSimThetaDegree = zeros(1,n);
AngleSimPhiDegree = zeros(1,n);

timeline = zeros(1, n);

for i=2:n
    %put a constant external force, ex: 0.3 pN
    fTemp = 0;
    if (i > appliedRange(1)*n && i < appliedRange(2)*n)
        fTemp = fLaser/10e9;
    else
        fTemp = 0;
    end
    timeline(i) = timeline(i-1) + deltat;
    extensionX = Xsim(i-1);
    extensionY = Ysim(i-1);
    extensionZ =(Zsim(i-1));
    
    anglePhi = AnglePhiSim(i-1);
    angleTheta = AngleThetaSim(i-1);

    extension = sqrt(extensionX^2 + extensionY^2 + extensionZ^2);
    FullExtension(i) = extension;
    DNAForceX = Marko_Sigga(kbT, Lp, Lo, extension, extensionX,'x', angleTheta, anglePhi);
    DNAForceY = Marko_Sigga(kbT, Lp, Lo, extension, extensionY,'y', angleTheta, anglePhi);
    DNAForceZ = Marko_Sigga(kbT, Lp, Lo, extension, extensionZ,'z', angleTheta, anglePhi);
    DNAForce(i) = 10e9 * sqrt(DNAForceX^2 + DNAForceY^2 + DNAForceZ^2);
    
    forceX = DNAForceX + fTemp*sin(alpha)*cos(beta); % giga Newton
    forceY = DNAForceY + fTemp*sin(alpha)*sin(beta); % giga  Newton
    forceZ = DNAForceZ + fTemp*cos(alpha); % giga  Newton
    
    
    
    %ldiff is change when tether is under confinement
    
    ldiff   = sqrt(2*D*deltat); 
    brownianTerms = normrnd(0,ldiff,3,1);
    brownianTermsX = brownianTerms(1);
    brownianTermsY = brownianTerms(2);
    brownianTermsZ = brownianTerms(3);
    
    sitaX = getDragCoefXY(Rb, extensionZ, sita0);
    sitaY = sitaX;
    sitaZ = getDragCoefZ(Rb, extensionZ, sita0);
        
    deltaX= (forceX*deltat)/sitaX + brownianTermsX;
    Xsim(i)=Xsim(i-1)+deltaX;%new element value  
   
    
    deltaY= (forceY*deltat)/sitaY + brownianTermsY;
    Ysim(i)=Ysim(i-1)+deltaY;%new element value 
    
    deltaZ= (forceZ*deltat)/sitaZ + brownianTermsZ;
    %if (Zsim(i-1) + deltaZ < Rb)
      %  Zsim(i)=Rb;%new element value 
   % else
        % Zsim(i)=Zsim(i-1)+deltaZ;%new element value 
    %end
    while (Zsim(i-1) + deltaZ < Rb)
        deltaZ = (forceZ*deltat)/sitaZ +  normrnd(0,ldiff,1,1);
    end
    Zsim(i)=Zsim(i-1)+deltaZ;
   
    
    
     AnglePhiSim(i) = atan(Ysim(i)/Xsim(i));
     AngleSimPhiDegree(i) = 180*(AnglePhiSim(i)/pi);
     
     AngleThetaSim(i) = atan(sqrt(Xsim(i)^2 + Ysim(i)^2) / Zsim(i));
     AngleSimThetaDegree(i) = 180*(AngleThetaSim(i)/pi);
    sita = 0;
    if (sitaX < sitaZ)
        sita = sitaX;
    else sita = sitaZ;
    end
     deltat = getNextTimestep(sita, Lp, Lo, extension, kbT);
  
    
end
