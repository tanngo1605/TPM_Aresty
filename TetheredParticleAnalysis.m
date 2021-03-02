function TetheredParticleAnalysis (Xdata,Ncorr,Nbins,deltat)
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
[Xsim, ForceSim, avgForce] = RandWalkSim(length(Xdata),deltat);  %generate simulated data
[Xsim_Original, ForceSim_Original, avgForce_Original] = RandWalkSim_Original(length(Xdata), deltat);
[FSim,rSim,histoSim] = GaussHistoX(Xsim,Nbins);
[FSim_Original,rSim_Original,histoSim_Original] = GaussHistoX(Xsim_Original,Nbins);
logACSim = LogAutoCorr (Xsim,Ncorr,deltat);
logACSim_Original = LogAutoCorr (Xsim_Original,Ncorr,deltat);

%%%%output%%%%
subplot(3,3,1);plot(time,Xdata,'b')%Experimental data
title ('Experimental data');
xlabel ('Time (sec)');
ylabel ('Bead Position (nm)');
subplot(3,3,6);
plot(time,ForceSim,'r')%Experimental data
str = sprintf('Force Over time (Average: %f)',avgForce);
title (str);
xlabel ('Time (sec)');
ylabel ('Force (Newton)');
subplot(3,3,2);plot(rData,FData,'bx',histoData(1,:),histoData(2,:),'bo')%Gaussian curve
title ('Gaussian fit (x) and Histogram (o) of exp data');
xlabel ('position (nm)');
ylabel ('Probability (1/nm)');
subplot(3,3,4);plot(time,Xsim,'r')%Simulated data
title ('Simulated data');
xlabel ('Time (sec)');
ylabel ('Bead Position (nm)');
subplot(3,3,5);plot(rSim,FSim,'rx',histoSim(1,:),histoSim(2,:),'ro')%Gaussian curve
title ('Gaussian fit (x) and Histogram (o)of Sim data');
xlabel ('position (nm)');
ylabel ('Probability (1/nm)');
subplot(3,3,7);plot(time,Xsim_Original,'g')%Simulated data _Original
title ('Simulated data with original model');
xlabel ('Time (sec)');
ylabel ('Bead Position (nm)');
subplot(3,3,8);plot(rSim_Original,FSim_Original,'gx',histoSim_Original(1,:),histoSim_Original(2,:),'go')%Gaussian curve_Original
title ('Gaussian fit (x) and Histogram (o)of Sim data with original model');
xlabel ('position (nm)');
ylabel ('Probability (1/nm)');
subplot(3,3,9);
plot(time,ForceSim_Original,'g')%Experimental data
str = sprintf('Force Over time (Average: %f)',avgForce_Original);
title (str);
xlabel ('Time (sec)');
ylabel ('Force (Newton)');
subplot(3,3,3);plot(logACData(1,:),logACData(2,:),'b-',logACSim(1,:),logACSim(2,:),'r-', logACSim_Original(1,:), logACSim_Original(2, :),'g-');%plot of correlation function of data
title ('Compare AutoCorrelation');
xlabel ('time difference (sec)');
ylabel ('Log[AutoCorrelation], (nm2)');

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

function force = Marko_Sigga (kbT, Lp, Lo, extension)
x = extension;
extension = abs(extension);
coff = (kbT)/(Lp);
% J /nm = J/ 10-9 m = J*10(9)/m 
change = 1/(4*(1 - extension/Lo)^2) - 1/4 + extension/Lo;
if (x>0)
    force = -coff*change;
else
    force = coff*change;
end

function newPos = validPosition(position)
    if (position > 1000) 
        newPos = 999;
    end
    if (position < -1000) 
        newPos= -999;
    end


function [Xsim, Fsim, avgForce]=RandWalkSim(n,deltat)
%This function simulates a 1D random walk in a harmonic potential
%%%%physical parameters%%%%
Lo       = 3477*.34;         % Lo = tether length (nm)
Lp      = 72;               % Lp = tether persistence length (nm)
Rb      = 240;              %bead radius (nm)
kbT     = 4.1*10^(-21);     %thermal energy (J)
eta     = 2.4*10^(-30);     %viscosity of H2O (J*s/nm^3)
D       = kbT/(6*pi*eta*Rb);%Stokes diffusion constant (nm2/s)
kappa   = 3/2*kbT/Lp/Lo;     %spring constant (J/nm2)
mu      = D*kappa/kbT*deltat;%bias of step in harmonic potential
ldiff   = sqrt(2*D*deltat); %diffusion length (nm)
deltaX = 0;
Xsim=zeros(1,n);
Fsim=zeros(1,n);
totalForce = 0;
for i=2:n
    %
    %k(x(n-1)) = F (Marko-Sigga)= f(x(n-1)) => k => mu => deltaX
    %extension = Lo + X(i)
    %forceMarko_Sigga = f (extension)
    deltaX= -Xsim(i-1)*mu + randn*ldiff;%step size for each element called
    Xsim(i)=Xsim(i-1)+deltaX;%new element value
    if ( abs(Xsim(i)) > 1000)
       %Xsim(i) = validPosition(Xsim(i));
       %Xsim
      
    end
       
        
    %Kappa = F_Marko_Sigga/deltax
   % kappa = Marko_Sigga(kbT, Lp, Lo, deltaX)/deltaX;
   % mu = D*kappa/kbT*deltat;
   extension = 0 + Xsim(i);
   
   %extension = Lp + deltaX;
   %extension = Lo + deltaX;
   %be careful with unit
   %create graph with MSForce over time
   %put a constant external force, ex: 0.3 pN
   
   force = Marko_Sigga(kbT, Lp, Lo, extension); %nN
   %%This is just the force applied along the DNA, not the net force.
   Fsim(i) = force*10e9; %Newton % F net
   %Xn = Xn+1 + (Fext + FBr)/y
   if (extension >=0)
   totalForce = totalForce + Fsim(i);
   else
   totalForce = totalForce - Fsim(i);
   end
   
   %kx + Fext = ma => x = (ma - Fext)/k
   kappa = -force/extension;
   mu = D*kappa/kbT*deltat;
    
end
avgForce = totalForce*10e12/n;

function [Xsim, Fsim, avgForce]=RandWalkSim_Original(n,deltat)
%This function simulates a 1D random walk in a harmonic potential
%%%%physical parameters%%%%
L       = 3477*.34;  %(1182.18) / 747; 1182.18 + 747       %tether length (nm)
Xi      = 72;               %tether persistence length (nm)
Rb      = 240;              %bead radius (nm)
kbT     = 4.1*10^(-21);     %thermal energy (J)
eta     = 2.4*10^(-30);     %viscosity of H2O (J*s/nm^3)
D       = kbT/(6*pi*eta*Rb);%Stokes diffusion constant (nm2/s)
kappa   = 3/2*kbT/Xi/L;     %spring constant (J/nm2)
mu      = D*kappa/kbT*deltat;%bias of step in harmonic potential
ldiff   = sqrt(2*D*deltat); %diffusion length (nm)
totalForce = 0;
Xsim=zeros(1,n);
Fsim=zeros(1,n);
for i=2:n
    deltaX=-Xsim(i-1)*mu + randn*ldiff;%step size for each element called
    Xsim(i)=Xsim(i-1)+deltaX;%new element value
    Fsim(i) = kappa*Xsim(i)*10e9;
    totalForce = totalForce + Fsim(i);
end
avgForce = totalForce*10e12/n;