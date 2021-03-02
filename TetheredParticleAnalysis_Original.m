function TetheredParticleAnalysis_Original (Xdata,Ncorr,Nbins,deltat)
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
Xdata     = transpose(Xdata); 
n         = length(Xdata);    %number of data points
time      = (1:n)*deltat;     %time series
[FData,rData,histoData] = GaussHistoX(Xdata,Nbins);
logACData = LogAutoCorr(Xdata,Ncorr,deltat);

%%%%Simulated data%%%%
Xsim = RandWalkSim(length(Xdata),deltat);  %generate simulated data
[FSim,rSim,histoSim] = GaussHistoX(Xsim,Nbins);
logACSim = LogAutoCorr (Xsim,Ncorr,deltat);

%%%%output%%%%
subplot(2,3,1);plot(time,Xdata,'b')%Experimental data
title ('Experimental data');
xlabel ('Time (sec)');
ylabel ('Bead Position (nm)');
subplot(2,3,2);plot(rData,FData,'bx',histoData(1,:),histoData(2,:),'bo')%Gaussian curve
title ('Gaussian fit (x) and Histogram (o) of exp data');
xlabel ('position (nm)');
ylabel ('Probability (1/nm)');
subplot(2,3,4);plot(time,Xsim,'r')%Simulated data
title ('Simulated data');
xlabel ('Time (sec)');
ylabel ('Bead Position (nm)');
subplot(2,3,5);plot(rSim,FSim,'rx',histoSim(1,:),histoSim(2,:),'ro')%Gaussian curve
title ('Gaussian fit (x) and Histogram (o)of Sim data');
xlabel ('position (nm)');
ylabel ('Probability (1/nm)');
subplot(2,3,3);plot(logACData(1,:),logACData(2,:),'b-',logACSim(1,:),logACSim(2,:),'r-')%plot of correlation function of data
title ('Compare AutoCorrelation');
xlabel ('time difference (sec)')
ylabel ('Log[AutoCorrelation], (nm2)')

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
    Xhisto(2,which)=Xhisto(2,which)+1;      %increment bin
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

function Xsim=RandWalkSim(n,deltat)
%This function simulates a 1D random walk in a harmonic potential
%%%%physical parameters%%%%
L       = 3477*.34;         %tether length (nm)
Xi      = 72;               %tether persistence length (nm)
Rb      = 240;              %bead radius (nm)
kbT     = 4.1*10^(-21);     %thermal energy (J)
eta     = 2.4*10^(-30);     %viscosity of H2O (J*s/nm^3)
D       = kbT/(6*pi*eta*Rb);%Stokes diffusion constant (nm2/s)
kappa   = 3/2*kbT/Xi/L;     %spring constant (J/nm2)
mu      = D*kappa/kbT*deltat;%bias of step in harmonic potential
ldiff   = sqrt(2*D*deltat); %diffusion length (nm)

Xsim=zeros(1,n);
for i=2:n
    deltaX=randn(1)*ldiff-Xsim(i-1)*mu;%step size for each element called
    Xsim(i)=Xsim(i-1)+deltaX;%new element value
end