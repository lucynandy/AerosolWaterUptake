function [aw] = WaterActivity(kappa,Ddry,RH)

    gf = GrowthFactor(kappa,Ddry,RH);
    
    aw = (g.^3 - 1)./(g.^3 - 1 + kappa);
    
end

function [gf] = GrowthFactor(kappa,Ddry,RH)

dg = 0.0000000001; % step size for derivative of function f
g1 = 1; % Give an initial starting growth factor guess.

% Newton raphson method
f = @(g) kappaCalc(Ddry,RH,g) - kappa;
df = @(g)(f(g+dg)-f(g))/dg;

for i = 1:10
    g2 = abs(g1 - f(g1)/df(g1));
    g1 = g2;
end

gf = g1;

end

function [kappa] = kappaCalc(Ddry,RH,g)

STw = 0.072;
Vw = 18*1e-6;
R = 8.314;
T = 298.15;

kappa = (g.^3 - 1).*(((exp((4.*STw.*Vw)./(R.*T.*g.*Ddry)))./(RH./100)) - 1);

end

