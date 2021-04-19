function [aw,gf,kappa] = WaterActivity(kappa1,Ddry,RH,STw,Vw,R,T)

[gf,kappa] = GrowthFactor(kappa1,Ddry,RH,STw,Vw,R,T);

aw = (gf.^3 - 1)./(gf.^3 - 1 + kappa1);

end

function [gf,kappa] = GrowthFactor(kappa1,Ddry,RH,STw,Vw,R,T)

dg = 0.0000000001; % step size for derivative of function f
g1 = 1; % Give an initial starting growth factor guess.

% Newton raphson method
f = @(g) kappaCalc(Ddry,RH,g,STw,Vw,R,T) - kappa1;
df = @(g)(f(g+dg)-f(g))/dg;

for i = 1:10
    g2 = abs(g1 - f(g1)/df(g1));
    g1 = g2;
end

gf = g1;

% % STw = 0.072; % N/m
% % Vw = 18e-6; % m^3/mol
% % R = 8.314; % J/mol-K
% % T = 298.15; % K
kappa = (gf.^3 - 1).*(((exp((4.*STw.*Vw)./(R.*T.*gf.*Ddry.*1e-6)))./(RH./100)) - 1);

end

function [kappa] = kappaCalc(Ddry,RH,g,STw,Vw,R,T)

% % STw = 0.072; % N/m
% % Vw = 18e-6; % m^3/mol
% % R = 8.314; % J/mol-K
% % T = 298.15; % K

kappa = (g.^3 - 1).*(((exp((4.*STw.*Vw)./(R.*T.*g.*Ddry.*1e-6)))./(RH./100)) - 1);

end

