function [aw,STsol,gf,RH,kappa] = WaterActivitySTsol(kappa_vf,Ddry,RH1,Vw,R,T)

aw = WaterActivity(kappa_vf,Ddry,RH1,Vw,R,T);

STsol = 0.001.*(159.25 - 85.146.*aw);
gf = (1 + ((kappa_vf.*aw)./(1 - aw))).^(1/3);
RH = 100.*aw.*exp((4.*STsol.*Vw)./(R.*T.*gf.*Ddry.*1e-6));
kappa = (gf.^3 - 1).*(((exp((4.*STsol.*Vw)./(R.*T.*gf.*Ddry.*1e-6)))./(RH./100)) - 1);

end

function aw = WaterActivity(kappa_vf,Ddry,RH1,Vw,R,T)

daw = 0.0000000001; % step size for derivative of function f
aw1 = 0.1; % Give an initial starting water activity guess.

% Newton raphson method
f = @(aw) RHcalc(aw,kappa_vf,Ddry,Vw,R,T) - RH1;
df = @(aw)(f(aw+daw)-f(aw))/daw;

for i = 1:10
    aw2 = abs(aw1 - f(aw1)/df(aw1));
    aw1 = aw2;
end

aw = aw1;

end

function RH = RHcalc(aw,kappa_vf,Ddry,Vw,R,T)

STsol = 0.001.*(159.25 - 85.146.*aw);
gf = (1 + ((kappa_vf.*aw)./(1 - aw))).^(1/3);
RH = 100.*aw.*exp((4.*STsol.*Vw)./(R.*T.*gf.*Ddry.*1e-6));

end

