function TetheredParticleAnalysis_WithDependentTimestep (Ncorr,Nbins, fLaser,phi, appliedRange,deltat, datapoints)
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
[Xsim, Ysim, FullExtension, AngleSimDegree] = RandWalkSim(n, fLaser,phi, appliedRange);  %generate simulated data
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

subplot(3,3,4);
plot(time, FullExtension,'r')%Experimental data
title ('Full extension');
xlabel ('Time (sec)');
ylabel ('Extension Bead Position (nm)');

subplot(3,3,5);
plot(time, AngleSimDegree,'b')%Experimental data
title ('Theta degree');
xlabel ('Time (sec)');
ylabel ('Degree (o)');

figure
subplot(1,1,1);
plot(Xsim, Ysim,'.')%Experimental data
title ('Position');
xlabel ('x');
ylabel ('y');


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



function force = Marko_Sigga (kbT, Lp, Lo, extension, direction, axis, angleTheta)
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
    tempForce = tempForce*abs(cos(angleTheta));
end
if (axis == 'y')
    tempForce = tempForce*abs(sin(angleTheta));
end
force = tempForce;

function newPos = validPosition(position)
    if (position > 1000) 
        newPos = 999;
    end
    if (position < -1000) 
        newPos= -999;
    end
    
function nextTimeStep = getNextTimestep(sita, Lp, Lo, x, kbT) %for 1 dim only
       constantCoeff = kbT/(Lp*Lo);
       change = 0.5/((1 - x/Lo)^3) + 1;
        gradient = abs(constantCoeff*change);
        nextTimeStep = (2*0.01 * sita) / gradient;

function [Xsim, Ysim, FullExtension, AngleSimDegree]=RandWalkSim(n,fLaser, phi,appliedRange)
%This function simulates a 1D random walk in a harmonic potential
%%%%physical parameters%%%%
Lo       = 3477*.34;            % Lo = tether length (nm) = 1182.18
Lp      = 72;                   % Lp = tether persistence length (nm)
Rb      = 240;                  %bead radius (nm)
kbT     = 4.1*10^(-21);         %thermal energy (J)
eta     = 2.4*10^(-30);         %viscosity of H2O (J*s/nm^3)
D       = kbT/(6*pi*eta*Rb);    %Stokes diffusion constant (nm2/s)
sita    = 6*pi*eta*Rb;          %drag coefficient, correct for close proximity for surface. FAXEN's law.
kappa   = 3/2*kbT/Lp/Lo;        %spring constant (J/nm2)
deltat = getNextTimestep(sita, Lp, Lo, 0, kbT);
Xsim=zeros(1,n);
Ysim=zeros(1,n);
FullExtension = zeros(1,n);
AngleSim = zeros(1,n);
AngleSimDegree = zeros(1,n);

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
    angleTheta = AngleSim(i-1);
    extension = sqrt(extensionX^2 + extensionY^2);
    FullExtension(i) = extension;
    
    forceX = Marko_Sigga(kbT, Lp, Lo, extension, extensionX,'x', angleTheta) + fTemp*cos(phi); % giga Newton
    forceY = Marko_Sigga(kbT, Lp, Lo, extension, extensionY,'y', angleTheta) + fTemp*sin(phi); % giga  Newton
    
    
    
    %ldiff is change when tether is under confinement
    
    ldiff   = sqrt(2*D*deltat); 
    brownianTerms = normrnd(0,ldiff,2,1);
    brownianTermsX = brownianTerms(1);
    brownianTermsY = brownianTerms(2);
    
    deltaX= (forceX*deltat)/sita + brownianTermsX;
    Xsim(i)=Xsim(i-1)+deltaX;%new element value  
   
    
    deltaY= (forceY*deltat)/sita + brownianTermsY;
    Ysim(i)=Ysim(i-1)+deltaY;%new element value 
    
     AngleSim(i) = atan(Ysim(i)/Xsim(i));
     %atan (Zsim(i) / (y^2 + x^2))
     AngleSimDegree(i) = 180*(AngleSim(i)/pi) ;
     deltat = getNextTimestep(sita, Lp, Lo, extension, kbT);
  
  
    
end
