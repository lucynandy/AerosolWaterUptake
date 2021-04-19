function [aw] = aw(m_total,molarRatio,Solute)

N = 2; % number of solutes

k = 1.38e-23;
T = 298.15;
Mw = 0.01802;  % Mol wt of H2O (kg/mol)
q1 = 1;
q2 = 2;
e = 1.60218e-19;
mu_w = 2.9;
r_ww = 2.82e-10;
D = 3.33564e-30;

for kk = 1:N
    if isequal(Solute(kk),{'GlutaricAcid'})
        j(kk) = 1; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 7; mu_j(kk) = 0.159751; r_jw(kk) = 1.03748E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
    elseif isequal(Solute(kk),{'SuccinicAcid'})
        j(kk) = 2; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 12; mu_j(kk) = 0.720837; r_jw(kk) = 1.72068E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
    elseif isequal(Solute,'CitricAcid')
        j(kk) = 3; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 9; mu_j(kk) = 16.1368; r_jw(kk) = 5E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
    elseif isequal(Solute(kk),{'sucrose'})
        j(kk) = 4; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 20; mu_j(kk) = 13.923; r_jw(kk) = 4.55659E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
    elseif isequal(Solute,'MaleicAcid')
        j(kk) = 5; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 3; mu_j(kk) = 4.3206; r_jw(kk) = 3.17995E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
    elseif isequal(Solute,'HNO3')
        j(kk) = 6; zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 37.287; n(kk) = 3; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 2.27E-10; r_pos(kk) = 2.1E-11; r_neg(kk) = 1.65E-10;
    elseif isequal(Solute,'NaNO3')
        j(kk) = 7;  zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 9.7517; n(kk) = 2; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 2.26E-10; r_pos(kk) = 1.01E-10; r_neg(kk) = 1.65E-10;
    elseif isequal(Solute,'NaOH')
        j(kk) = 8; zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 11.706; n(kk) = 4; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 1.5E-10; r_pos(kk) = 1.01E-10; r_neg(kk) = 1.53E-10;
    elseif isequal(Solute,'NaCl')
        j(kk) = 9; zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 8.5709; n(kk) = 4; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 2.24E-10; r_pos(kk) = 1.01E-10; r_neg(kk) = 1.81E-10;
    elseif isequal(Solute,'NH4NO3')
        j(kk) = 10;  zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 8.0672; n(kk) = 2; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 2.16E-10; r_pos(kk) = 1.78E-10; r_neg(kk) = 1.65E-10;
    elseif isequal(Solute,'NH4Cl')
        j(kk) = 11;  zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 14.586; n(kk) = 3; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 2.23E-10; r_pos(kk) = 1.78E-10; r_neg(kk) = 1.81E-10 ;
    elseif isequal(Solute,'HCl')
        j(kk) = 12;  zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 28.172; n(kk) = 7; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 1.47E-10; r_pos(kk) = 0.21E-10; r_neg(kk) = 1.81E-10 ;
    elseif isequal(Solute,'NaHCO3')
        j(kk) = 13;  zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 13.00; n(kk) = 2; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 4.02E-10; r_pos(kk) = 1.02E-10; r_neg(kk) = 1.56E-10 ;
    elseif isequal(Solute(kk),{'NH42SO4'})
        j(kk) = 14; zpos(kk) = 1; zneg(kk) = 2; v(kk) = 3; rho(kk) = 9.5964; n(kk) = 3; mu_j(kk) = 5.4394; r_jw(kk) = 3.4109e-10; r_jj(kk) = 0; r_pos(kk) = 1.37e-10; r_neg(kk) = 2.58e-10 ;
    elseif isequal(Solute,'Na2SO4')
        j(kk) = 15;  zpos(kk) = 1; zneg(kk) = 2; v(kk) = 3; rho(kk) = 8.110; n(kk) = 3; mu_j(kk) = 47.02; r_jw(kk) = 7.32e-10; r_jj(kk) = 4.3448e-10; r_pos(kk) = 1.02e-10; r_neg(kk) = 2.58e-10 ;
    elseif isequal(Solute,'H2CO3')
        j(kk) = 16;  zpos(kk) = 1; zneg(kk) = 2; v(kk) = 3; rho(kk) = 13.00; n(kk) = 2; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 4.02E-10; r_pos(kk) = 0.21E-10; r_neg(kk) = 1.79E-10 ;
    end
    
    for ii = 1:n(kk)-1
        deltaE(kk,ii) = ((mu_j(kk)*mu_w*D^2)/((1.113e-10)*(r_jw(kk)+(ii-1)*r_ww)^3))-((mu_w^2*D^2)/((1.113e-10)*(ii*r_ww)^3));
        C(kk,ii) = exp(deltaE(kk,ii)/(k*T)); % Energy parameter
    end
    
end

for kk = 1:N-1
    soluteratio(kk) = molarRatio;
end

assignin('base', 'soluteratio', soluteratio);

t = length(m_total);

for i = 1:t
    aw(N,i) = NR(j,rho,v,n,zpos,zneg,C,m_total,N,soluteratio);
end

end

function [awOUT] = NR(j,rho,v,n,zpos,zneg,C,m_total,N,soluteratio)

for i = 1:t
        m(N,i) = m_total(i)./(soluteratio.*v(N-1) + v(N));
        m(N-1,i) = soluteratio.*m(N,i);
end

daw = 0.0000000001; % step size for derivative of function f
aw1 = 1; % Give an initial starting molality guess.

% Newton raphson method
f = @(aw) awCalc(j,rho,v,n,zpos,zneg,C,m(N,:),aw,N,soluteratio) - aw;
df = @(aw)(f(aw+daw)-f(aw))/daw;
for i = 1:10
    aw2 = abs(aw1 - f(aw1)/df(aw1));
    aw1 = aw2;
end

awOUT = aw1;

end

function [OUT] = awCalc(j,rho,v,n,zpos,zneg,C,m,aw,N,soluteratio)

Mw=0.0180152; % Molecular weight water (kg/mol)
Ax=2.917; % Debye-Huckel coeff. mol frac @ 298.15K

% disp(aw); fprintf('\b'); disp(m); fprintf('\b'); disp(soluteratio);

IxTop = zpos(N)*zneg(N).*v(N);
for kk = 1:N-1
    IxTop = (IxTop + soluteratio(kk)*zpos(kk)*zneg(kk).*v(kk));
    IxTop = m*IxTop;
end

IxBottom = v(N);
for kk = 1:N-1
    IxBottom = (IxBottom + soluteratio(kk).*v(kk));
    IxBottom = m*IxBottom;
end

Ix = (1/2)* (IxTop)./(IxBottom+1/Mw);

KwTop = (zpos(N)*zneg(N)*v(N))./(1+rho(N)*Ix.^0.5);

for kk = 1:N-1
    KwTop = (KwTop + (zpos(kk)*zneg(kk).*v(kk)*soluteratio(kk))./(1+rho(kk)*Ix.^0.5));
    KwTop = m*KwTop;
end

Kw = exp((Ax.*(Ix.^0.5).*(KwTop))./(IxBottom+1/Mw));

awbar=aw./(Kw);

% Calculate molality Using EQN.27
for p =1:N
    
    if j(p) >= 1 && j(p) <= 5
        NumorSum(p)=0;
        for k=1:(n(p)-1)
            NumorSum(p)=NumorSum(p)+(awbar.^k).*(1-C(p,k)).*prod(C(p,1:k-1));
        end
        DenomSum(p)=0;
        for k=1:(n(p)-2)
            DenomSum(p)=DenomSum(p)+k*(awbar.^(k-1)).*prod(C(p,1:k));
        end
        Denom(p)=(1-awbar).^2.*(DenomSum(p))+(n(p)-1-(n(p)-2).*awbar).*awbar.^(n(p)-2).*prod(C(p,1:n(p)-1));
        m0(p)=((1-awbar)/(Mw*v(p)*awbar))*(1-NumorSum(p))/Denom(p);
    else
        NumorSum(p)=0;
        for k=1:(n(p)-1)
            NumorSum(p)=NumorSum(p)+(awbar.^k).*(1-C(p,k)).*prod(C(p,1:k-1));
        end
        DenomSum(p)=0;
        for k=1:(n(p)-2)
            DenomSum(p)=DenomSum(p)+k*(awbar.^(k-1)).*prod(C(p,1:k));
        end
        Denom(p)=(1-awbar).^2.*(DenomSum(p))+(n(p)-1-(n(p)-2).*awbar).*awbar.^(n(p)-2).*prod(C(p,1:n(p)-1));
        m0(p)=((1-awbar)/(Mw*v(p)*awbar))*(1-NumorSum(p))/Denom(p);
    end
    
end

mixturemodel=1./m0(N);

for k=1:N-1
    mixturemodel=mixturemodel+soluteratio(kk)/m0(kk);
end

m = 1./mixturemodel;

OUT = awbar.*Kw;

end

