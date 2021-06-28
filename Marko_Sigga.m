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