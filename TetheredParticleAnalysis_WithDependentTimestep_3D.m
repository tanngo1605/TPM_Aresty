%expand to 3dim
%change the Marko SIgga model

function TetheredParticleAnalysis_WithDependentTimestep_3D (Ncorr,Nbins, fLaser,alpha, beta,appliedRange,deltat,datapoints)
%theta: angle with z axis
%phi: angle in x0y plane
%alpha: angle of force with z axis
%beta: angle of the force in xOy plane.


%This function simulates and analyzes the motion of a 1D random walk
%confined in a harmonic potential well.
%
%Written by Luke Sullivan, Ursinus College
%Edited by John F. Beausang, University of Pennsylvania
%
%Xdata    = array of bead position (nm)
%nCorr    = number of points in correlation function
%Nbins    = number of histogram bins
%deltat   = time step of data and simulation (sec)

%%%%Experimental Data%%%%
%Xdata     = transpose(Xdata); 
n         = datapoints;    %number of data points
time      = (1:n)*deltat;     %time series
%[FData,rData,histoData] = GaussHistoX(Xdata,Nbins);
%logACData = LogAutoCorr(Xdata,Ncorr,deltat);

%%%%Simulated data%%%%
[Xsim, Ysim,Zsim, FullExtension, AngleThetaSimDegree, AnglePhiSimDegree] = RandWalkSim(n, fLaser,alpha, beta, appliedRange);  %generate simulated data
%[FSimX,rSimX,histoSimX] = GaussHistoX(Xsim,Nbins);
%[FSimY,rSimY,histoSimY] = GaussHistoX(Ysim,Nbins);

%logACSim = LogAutoCorr (Xsim,Ncorr,deltat);

%%%%output%%%%
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

figure
subplot(1,1,1);
plot3(Xsim, Ysim, Zsim,'.')%Experimental data
title ('Position');
xlabel ('x');
ylabel ('y');
zlabel('z')
grid on


%%%%Subroutines%%%%%

function [F,r,Xhisto]=GaussHistoX (Xdata,Nbins)
%This function histograms the data and fits a Gaussian distribution
Xmax        = max(abs(Xdata));  %maximum position
n           = length(Xdata);    %number of data points
binWidth    = Xmax/Nbins;       %histogram bin width
stdevX      = std(Xdata);       %standard deviation of data
F=zeros(1,n);r=zeros(1,n);Xhisto=zeros(2,Nbins+1);  %initialize
Xhisto(1,:)= ((1:Nbins+1)-.5)*binWidth;     %midpoint of histogram bins
for i=1:n
    r(i)=abs(Xdata(i));
    F(i)=2/sqrt(2*pi*stdevX^2)*exp(-(Xdata(i)^2)/(2*stdevX^2));%1 sided gaussian curve
    which=1+floor(abs(Xdata(i))/binWidth);  %which bin data falls into
    if (which > 0)
    Xhisto(2,which)=Xhisto(2,which)+1;%increment bin
    end
    %temporary ignore
end
Xhisto(2,:)=Xhisto(2,:)/n/binWidth; %convert counts to probability (1/nm)

function logac = LogAutoCorr (Xdata,Ncorr,deltat)
%This function determines the autocorrelation of the data for Ncorr points
n         = length(Xdata);
logac     = zeros(2,Ncorr);
logac(1,:)= (0:Ncorr-1)*deltat;      %time steps
for s = 1:Ncorr
    temp = zeros (1,n-s+1);
    for i=1:(n-s+1)
        temp(i)=Xdata(i)*Xdata(i+s-1);
    end
    logac(2,s)=log10(sum(temp)/(n-s+1));
end



function force = Marko_Sigga (kbT, Lp, Lo, extension, direction, axis, angleTheta, anglePhi)
extension = abs(extension);
coff = (kbT)/(Lp)*1.0;
% J /nm = J/ 10-9 m = J*10(9)/m 
something = 4.0*(1-extension/Lo)*(1-extension/Lo);
%change = 1/(4*(1 - extension/Lo)^2) - 1/4 + extension/Lo;
change = 1.0/something - 1.0/4 + 1.0*extension/Lo;
%Should be < 65pN
%tempForce = coff*change;
tempForce = 0;
if (direction>0)
    tempForce = -1.0*coff*change;
else
    tempForce = 1.0*coff*change;
end
if (axis == 'x')
    tempForce = tempForce*abs(sin(angleTheta)*cos(anglePhi));
end
if (axis == 'y')
    tempForce = tempForce*abs(sin(angleTheta)*sin(anglePhi));
end
if (axis == 'z')
    tempForce = tempForce*abs(cos(angleTheta));
end
force = tempForce;

function newPos = validPosition(position)
    if (position > 1000) 
        newPos = 999;
    end
    if (position < -1000) 
        newPos= -999;
    end
    
function nextTimeStep = getNextTimestep(sita, Lp, Lo, x, kbT) 
       constantCoeff = kbT/(Lp*Lo);
       change = 0.5/((1 - x/Lo)^3) + 1;
        gradient = abs(constantCoeff*change);
        nextTimeStep = (2*0.01 * sita) / gradient;

function [Xsim, Ysim, Zsim, FullExtension, AngleSimThetaDegree , AngleSimPhiDegree]=RandWalkSim(n,fLaser,alpha, beta,appliedRange)
%This function simulates a 1D random walk in a harmonic potential
%%%%physical parameters%%%%
Lo       = 3477*.34;            % Lo = tether length (nm) = 1182.18
Lp      = 72;                   % Lp = tether persistence length (nm)
Rb      = 240;                  %bead radius (nm)
kbT     = 4.1*10^(-21);         %thermal energy (J)
eta     = 2.4*10^(-30);         %viscosity of H2O (J*s/nm^3)
D       = kbT/(6*pi*eta*Rb);    %Stokes diffusion constant (nm2/s)
sita    = 6*pi*eta*Rb;          %drag coefficient, correct for close proximity for surface. FAXEN's law, 
kappa   = 3/2*kbT/Lp/Lo;        %spring constant (J/nm2)
deltat = getNextTimestep(sita, Lp, Lo, 0, kbT);
Xsim=zeros(1,n);
Ysim=zeros(1,n);
Zsim=zeros(1,n);

FullExtension = zeros(1,n);

AnglePhiSim = zeros(1,n);

AngleThetaSim = zeros(1,n);
AngleThetaSim = AngleThetaSim + pi/2;



AngleSimThetaDegree = zeros(1,n);


AngleSimPhiDegree = zeros(1,n);

for i=2:n
    %put a constant external force, ex: 0.3 pN
    fTemp = 0;
    if (i > appliedRange(1)*n && i < appliedRange(2)*n)
        fTemp = fLaser/10e9;
    else
        fTemp = 0;
    end
    extensionX = Xsim(i-1);
    extensionY = Ysim(i-1);
    extensionZ =(Zsim(i-1));
    
    anglePhi = AnglePhiSim(i-1);
    angleTheta = AngleThetaSim(i-1);
    
    extension = sqrt(extensionX^2 + extensionY^2 + extensionZ^2);
    FullExtension(i) = extension;
    
    forceX = Marko_Sigga(kbT, Lp, Lo, extension, extensionX,'x', angleTheta, anglePhi) + fTemp*sin(alpha)*cos(beta); % giga Newton
    forceY = Marko_Sigga(kbT, Lp, Lo, extension, extensionY,'y', angleTheta, anglePhi) + fTemp*sin(alpha)*sin(beta); % giga  Newton
    forceZ = Marko_Sigga(kbT, Lp, Lo, extension, extensionZ,'z', angleTheta, anglePhi) + fTemp*cos(alpha); % giga  Newton
    
    
    
    %ldiff is change when tether is under confinement
    
    ldiff   = sqrt(2*D*deltat); 
    brownianTerms = normrnd(0,ldiff,3,1);
    brownianTermsX = brownianTerms(1);
    brownianTermsY = brownianTerms(2);
    brownianTermsZ = brownianTerms(3);
    
    deltaX= (forceX*deltat)/sita + brownianTermsX;
    Xsim(i)=Xsim(i-1)+deltaX;%new element value  
   
    
    deltaY= (forceY*deltat)/sita + brownianTermsY;
    Ysim(i)=Ysim(i-1)+deltaY;%new element value 
    
    deltaZ= (forceZ*deltat)/sita + brownianTermsZ;
    Zsim(i)=Zsim(i-1)+deltaZ;%new element value 
    
    
     AnglePhiSim(i) = atan(Ysim(i)/Xsim(i));
     AngleSimPhiDegree(i) = 180*(AnglePhiSim(i)/pi);
     
     AngleThetaSim(i) = atan(sqrt(Xsim(i)^2 + Ysim(i)^2) / Zsim(i));
     AngleSimThetaDegree(i) = 180*(AngleThetaSim(i)/pi);
   
     deltat = getNextTimestep(sita, Lp, Lo, extension, kbT);
  
  
    
end
