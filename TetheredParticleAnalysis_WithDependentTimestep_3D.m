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
[Xsim, Ysim,Zsim, FullExtension, AngleThetaSimDegree, AnglePhiSimDegree, DNAForce] = RandWalkSim(n, fLaser,alpha, beta, appliedRange);  %generate simulated data
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

subplot(3,3,7);
plot(time, DNAForce,'b')%Experimental data
title ('Stretching DNA force');
xlabel ('Time (sec)');
ylabel ('Magnitude (N)');

% figure
% subplot(1,1,1);
% for timeStep = 1:n
%    if (timeStep > appliedRange(1)*n && timeStep < appliedRange(2)*n)
%        plot3(Xsim(timeStep), Ysim(timeStep), Zsim(timeStep), 'r.');
%    else
%        plot3(Xsim(timeStep), Ysim(timeStep), Zsim(timeStep), 'b.');
%    end
%    pause(0.05);
%    hold on;
% end
% title ('Position');
% xlabel ('x');
% ylabel ('y');
% zlabel('z')
% grid on


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
d=1.6;
kbT = kbT * 10^21;
Ko = 16*kbT * Lp * d^-2;
extension = abs(extension);
func=@(extension, Fe) kbT/Lp*(1/4*1/(1-extension/Lo+Fe/Ko)^2-1/4+extension/Lo-Fe/Ko)-Fe;
tempFunc = @(Fe) func(extension, Fe);
tempForce = fzero(tempFunc,63)/10^21;

if (direction>0)
    tempForce = -1.0*tempForce;
else
    tempForce = 1.0*tempForce;
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


% function force = Marko_Sigga (kbT, Lp, Lo, extension, direction, axis, angleTheta, anglePhi)
% extension = abs(extension);
% coff = (kbT)/(Lp)*1.0;
% % J /nm = J/ 10-9 m = J*10(9)/m 
% something = 4.0*(1-extension/Lo)*(1-extension/Lo);
% %change = 1/(4*(1 - extension/Lo)^2) - 1/4 + extension/Lo;
% change = 1.0/something - 1.0/4 + 1.0*extension/Lo;
% %Should be < 65pN
% %tempForce = coff*change;
% tempForce = 0;
% if (direction>0)
%     tempForce = -1.0*coff*change;
% else
%     tempForce = 1.0*coff*change;
% end
% if (axis == 'x')
%     tempForce = tempForce*abs(sin(angleTheta)*cos(anglePhi));
% end
% if (axis == 'y')
%     tempForce = tempForce*abs(sin(angleTheta)*sin(anglePhi));
% end
% if (axis == 'z')
%     tempForce = tempForce*abs(cos(angleTheta));
% end
% force = tempForce;
% 
% function newPos = validPosition(position)
%     if (pos/ition > 1000) 
%         newPos = 999;
%     end
%     if (position < -1000) 
%         newPos= -999;
%     end
    
function nextTimeStep = getNextTimestep(sita, Lp, Lo, x, kbT) 
       constantCoeff = kbT/(Lp*Lo);
       change = 0.5/((1 - x/Lo)^3) + 1;
        gradient = abs(constantCoeff*change);
        nextTimeStep = (2*0.01 * sita) / gradient;
        
function dragCoefXY = getDragCoefXY(r, z, sita0)
ratio = r/z;
denom = 1 - (9/16)*ratio + (1/8)*ratio^3 - (45/256)*ratio^4 - (1/16)*ratio^5;
dragCoefXY = sita0 / denom;
function dragCoefZ = getDragCoefZ(r, z, sita0)
ratio = r/z;
denom = 1 - (9/8)*ratio + 0.5*ratio^3 - 0.57*ratio^4 + 0.2*ratio^5 + (7/200)*ratio^11 - (1/25)*ratio^12;
dragCoefZ = sita0 / denom;
% function nextTimeStep = getNextTimestep(kbT, sita, Lp, Lo, extension, extensionX, extensionY, extensionZ, curForce )
%     kbT = kbT * 10^21;
%     d=1.6;
%     Ko = 16*kbT * Lp * d^-2;
%     myFunc = @(extension, componentExtension,curForce, devF) kbT/Lp * ( -0.5* 1/(1 - extension/Lo + curForce/Ko )^3 * (1 - componentExtension/(extension*Lo) + devF/Ko) + componentExtension/(extension*Lo) - devF/Ko) - devF;
%     
%     tempFunc = @(devF) myFunc(extension, extensionX, curForce, devF);
%     devF_X = fzero(tempFunc,0);
%     
%     tempFunc = @(devF) myFunc(extension, extensionY, curForce, devF);
%     devF_Y = fzero(tempFunc,0);
%     
%     tempFunc = @(devF) myFunc(extension, extensionZ, curForce, devF);
%     devF_Z = fzero(tempFunc,0);
%     
%     gradient = sqrt(devF_X^2 + devF_Y^2 + devF_Z^2);
%     
%     
%     nextTimeStep = (2*0.01 * sita) / gradient;

function [Xsim, Ysim, Zsim, FullExtension, AngleSimThetaDegree , AngleSimPhiDegree, DNAForce]=RandWalkSim(n,fLaser,alpha, beta,appliedRange)
%This function simulates a 1D random walk in a harmonic potential
%%%%physical parameters%%%%
Lo       = 3477*.34;            % Lo = tether length (nm) = 1182.18
Lp      = 72;                   % Lp = tether persistence length (nm)
Rb      = 240;                  %bead radius (nm)
Rb = 50;
%change to dynamic kbT = k * T
kbT     = 4.1*10^(-21);         %thermal energy (J)
eta     = 2.4*10^(-30);         %viscosity of H2O (J*s/nm^3)
D       = kbT/(6*pi*eta*Rb);    %Stokes diffusion constant (nm2/s)
sita0    = 6*pi*eta*Rb;          %drag coefficient, correct for close proximity for surface. FAXEN's law, 
kappa   = 3/2*kbT/Lp/Lo;        %spring constant (J/nm2)
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
%      if (extension > 0)
%      deltat = getNextTimestep(kbT, sita, Lp, Lo, extension, abs(extensionX),abs(extensionY),abs(extensionZ), sqrt(DNAForceX^2 + DNAForceY^2 + DNAForceZ^2)*10^21);
%      else
%          deltat = 0.0001;
%      end
%      time(i) = time(i-1) + deltat;
  
  
    
end
