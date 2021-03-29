clear all;
close all;
kbT = 4.1*10^-21;
Lo       = 3477*.34;            % Lo = tether length (nm) = 1182.18
Lp      = 72;                   % Lp = tether persistence length (nm);
d=2;
%Ko = 16*kbT * Lp/(d^2)*10^21; %GN 
Ko = 16*kbT * Lp/(d^2);
%Ko= 1100;
% extensions = 0:10:1180;
% 
% force1 = zeros(1, length(extensions));
% force2 = zeros(1, length(extensions));
% 
% for i = 1: length(extensions)
%     extension = extensions(i);
%     syms F positive
%     eqn1 = F == (kbT/Lp) * ( 0.25*(1 - extension/Lo)^-2 - 0.25 + extension/Lo);
%     eqn2 = 0 == (kbT/Lp) * ( 0.25*(1 - extension/Lo + F/Ko)^-2 - 0.25 + extension/Lo - F/Ko) - F;
%     %J/nm -> J/ -> 10^9 N 
%     % J/m -> N
%    % force1(i) = vpasolve(eqn1, F)*10^21;
% %assume(F,'real');
%     force2(i) = vpasolve(eqn2, F, [0 Inf])*10^21;
% end
% plot(extensions, force1, 'r');
% hold on
% plot(extensions, force2, 'b');
% axis([0, 1300, 0, 60]);

% %syms F
% syms F positive
% eqn1 = F == (kbT/Lp) * ( 0.25*(1 - extension/Lo)^-2 - 0.25 + extension/Lo);
% eqn2 = F == (kbT/Lp) * ( 0.25*(1 - extension/Lo + F/Ko)^-2 - 0.25 + extension/Lo - F/Ko);
% 
% format long
% %force = solve(eqn, F, 'Real', true,'IgnoreAnalyticConstraints',true)*1.0;
% force = vpasolve(eqn1, F);
% %assume(F,'real');
% force2 = vpasolve(eqn2, F, [0 Inf]);
%func=@(extension, F)((kbT/Lp) * ( 0.25.*(1 - extension./Lo + F./Ko).^-2 - 0.25 + extension./Lo - F./Ko) - F);
%func=@(extension, F) extension.^2 + F.^2 - 1;

%% Working Plot
% func=@(extension, Fe) kbT/Lp*(1/4*1/(1-extension./Lo+Fe./Ko)^2-1/4+extension./Lo-Fe./Ko)-Fe;
% figure
%     hold
%     fimplicit(func, [0 1300 0 70*10^-21])
%     plot([Lo Lo],[0 100*10^-21],'k')
%     xlabel('Extension [-]')
%     ylabel('F [GN]')
%     title('WLC modified Marko & Siggia  1995')

%% try vpasolve
% % extension = Lo * 1.04
% %extensions = 0:10:(Lo*1.05);
% extensions = 0:1:1230;
% force = zeros(1, length(extensions));
% %myfun = @(extension) kbT/Lp*(1/4*1/(1-extension/Lo+Fe/Ko)^2-1/4+extension/Lo-Fe/Ko)-Fe;  % parameterized function
% tenp = 0;
% for i= 1 : length(extensions)
%     syms Fe
%     extension = extensions(i);
%     %fun = myfun(extension);
%     eqn = kbT/Lp*(1/4*1/(1-extension/Lo+Fe/Ko)^2-1/4+extension/Lo-Fe/Ko)-Fe == 0; 
%     force(i) = vpasolve(kbT/Lp*(1/4*1/(1-extension/Lo+Fe/Ko)^2-1/4+extension/Lo-Fe/Ko)-Fe == 0, Fe, [-Inf Inf]);
%     temp = force(i);
%     %force(i) = fzero(fun, 65*10^-21);
%     
% end
% 
% plot(extensions, force);
% hold on
% plot([Lo, Lo], [0, 70*10^-21]);

%% try fzero, uplift the eqn by 10^21, find in interval [0, 65]
Lp = 72; %Contour length [nm]
kbT=4.1; %[pN*nm]
d=2;
Ko = (16*kbT * Lp/(d^2));

extensions = 0:1:1250;
func=@(extension, Fe) kbT/Lp*(1/4*1/(1-extension/Lo+Fe/Ko)^2-1/4+extension/Lo-Fe/Ko)-Fe;
force = zeros(1, length(extensions));
for i = 1 : length(extensions)
    tempFunc = @(Fe) func(extensions(i), Fe);
    force(i) = fzero(tempFunc, 0); %pN
end

plot(extensions, force);
hold on;
plot([Lo, Lo], [0, 65], 'r');
