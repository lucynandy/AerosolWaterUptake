function [ m, Ix, a, f_Inf, f_FS, gamma_Inf ] = MolalityMixture(aw,molarRatio,Solute,N)

% N = 2; % number of solutes

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
    
%     prompt = 'Enter solute of aqueous solution     ';
%     Solute = input(prompt,'s');
  
    if isequal(Solute(kk),{'GlutaricAcid'})
        j(kk) = 1; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 7; mu_j(kk) = 0.159751; r_jw(kk) = 1.03748E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;

%         dataloc = 'C:\Postdoctoral research\Dicarboxylic Acids\DissociationModel\Data\';
%         datafile = strcat(dataloc,Solute,'.m');
%         run(datafile)

% % % %         m_data(:,kk) = x     % x = molality
% % % %         o_data(:,kk) = Q    % Q = Osmotic coefficient
% % % %         aw_data(:,kk) = exp(-o_data(:,kk).*v.*Mw.*m_data(:,kk))
% % % %         x_data(:,kk) = (v.*m_data(:,kk))./(v.*m_data(:,kk)+1/Mw)

%         m_data(:,kk) = zeros(size(x));
%         o_data(:,kk) = zeros(size(x));
%         aw_data(:,kk) = zeros(size(x));
%         x_data(:,kk) = zeros(size(x));
        
%         m_data(kk,:) = x';
%         o_data(kk,:) = Q';
%         aw_data(kk,:) = exp(-Mw*v.*Q.*x);
%         x_data(kk,:) = (v.*x)./(v.*x+1/Mw);

%         if isequal(Solute,'MalicAcid')
        %         j(kk) = 1; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 4; mu_j(kk) = 7.67324; r_jw(kk) = 3.57E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
        %         if isequal(Solute,'ButyricAcid')
        %         j(kk) = 1; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 3; mu_j(kk) = 157.781; r_jw(kk) = 13.5E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
        %         if isequal(Solute,'AceticAcid')
        %         j(kk) = 1; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 6; mu_j(kk) = 1.43626; r_jw(kk) = 2.28E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
    elseif isequal(Solute(kk),{'SuccinicAcid'})
        j(kk) = 2; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 12; mu_j(kk) = 0.720837; r_jw(kk) = 1.72068E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
%         mu_j(kk) = 0.1208*(r_jw(kk)*10^10)^3;
%         mu_j(kk) = 2.70024; r_jw(kk) = 2.81706E-10;
        
%         dataloc = 'C:\Postdoctoral research\Dicarboxylic Acids\DissociationModel\Data\';
%         datafile = strcat(dataloc,Solute,'.m');
%         run(datafile)

%         m_data = x;     % x = molality
%         o_data = Q;    % Q = Osmotic coefficient
%         aw_data = exp(-o_data.*v.*Mw.*m_data);
%         x_data = (v.*m_data)./(v.*m_data+1/Mw);
% 

        % elseif isequal(Solute,'TartaricAcid')
        %                     j(kk) = 2; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 4; mu_j(kk) = 77.68; r_jw(kk) = 8.63E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
%             elseif isequal(Solute,'glycerol')
%                 j(kk) = 2; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 12; mu_j(kk) = 3.804; r_jw(kk) = 3.04E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
    elseif isequal(Solute,'CitricAcid')
        j(kk) = 3; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 9; mu_j(kk) = 16.1368; r_jw(kk) = 5E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
    elseif isequal(Solute(kk),{'sucrose'})
        j(kk) = 4; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 20; mu_j(kk) = 13.923; r_jw(kk) = 4.55659E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
    elseif isequal(Solute,'MaleicAcid')
        j(kk) = 5; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 3; mu_j(kk) = 4.3206; r_jw(kk) = 3.17995E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
%     elseif isequal(Solute,'MalonicAcid')
%         j(kk) = 5; zpos(kk) = 0; zneg(kk) = 0; v(kk) = 1; rho(kk) = 0; n(kk) = 3; mu_j(kk) = 40.2079; r_jw(kk) = 7.09E-10; r_jj(kk) = 0; r_pos(kk) = 0; r_neg(kk) = 0;
    elseif isequal(Solute,'HNO3')
        j(kk) = 6; zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 37.287; n(kk) = 3; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 2.27E-10; r_pos(kk) = 2.1E-11; r_neg(kk) = 1.65E-10;
    elseif isequal(Solute,'NaNO3')
        j(kk) = 7;  zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 9.7517; n(kk) = 2; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 2.26E-10; r_pos(kk) = 1.01E-10; r_neg(kk) = 1.65E-10;
        %     elseif isequal(Solute,'LiCl')
        %         j(kk) = 7; zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 27.275; n(kk) = 6; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 1.25E-10; r_pos(kk) = 7.6E-11; r_neg(kk) = 1.81E-10;
    elseif isequal(Solute,'NaOH')
        j(kk) = 8; zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 11.706; n(kk) = 4; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 1.5E-10; r_pos(kk) = 1.01E-10; r_neg(kk) = 1.53E-10;
    elseif isequal(Solute(kk),{'NaCl'})
        j(kk) = 9; zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 8.5709; n(kk) = 4; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 2.24E-10; r_pos(kk) = 1.01E-10; r_neg(kk) = 1.81E-10;
%         dataloc = 'C:\Postdoctoral research\DataLit\Data_Single_Salts\';
%         datafile = strcat(dataloc,Solute,'.m');
%         run(datafile)

% % % %         m_data = x;     % x = molality
% % % %         o_data = Q;    % Q = Osmotic coefficient
% % % %         aw_data = exp(-o_data.*v.*Mw.*m_data);
% % % %         x_data = (v.*m_data)./(v.*m_data+1/Mw);

%         m_data(:,kk) = zeros(size(x));
%         o_data(:,kk) = zeros(size(x));
%         aw_data(:,kk) = zeros(size(x));
%         x_data(:,kk) = zeros(size(x));
        
%         m_data(kk) = x';
%         o_data(kk) = Q';
%         aw_data(kk) = exp(-Mw*v.*Q.*x);
%         x_data(kk) = (v.*x)./(v.*x+1/Mw);
        
    elseif isequal(Solute,'NH4NO3')
        j(kk) = 10;  zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 8.0672; n(kk) = 2; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 2.16E-10; r_pos(kk) = 1.78E-10; r_neg(kk) = 1.65E-10;
    elseif isequal(Solute,'NH4Cl')
        j(kk) = 11;  zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 14.586; n(kk) = 3; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 2.23E-10; r_pos(kk) = 1.78E-10; r_neg(kk) = 1.81E-10 ;
    elseif isequal(Solute,'HCl')
        j(kk) = 12;  zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 28.172; n(kk) = 7; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 1.47E-10; r_pos(kk) = 0.21E-10; r_neg(kk) = 1.81E-10 ;
    elseif isequal(Solute,'NaHCO3')
        j(kk) = 13;  zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 13.00; n(kk) = 2; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 4.02E-10; r_pos(kk) = 1.02E-10; r_neg(kk) = 1.56E-10 ;
%     elseif isequal(Solute,'Na2CO3')
%         j(kk) = 14;  zpos(kk) = 1; zneg(kk) = 2; v(kk) = 3; rho(kk) = 13.60; n(kk) = 2; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 4.02E-10; r_pos(kk) = 1.02E-10; r_neg(kk) = 1.79E-10 ;
    elseif isequal(Solute(kk),{'NH42SO4'})
        j(kk) = 14; zpos(kk) = 1; zneg(kk) = 2; v(kk) = 3; rho(kk) = 9.5964; n(kk) = 3; mu_j(kk) = 5.4394; r_jw(kk) = 3.4109e-10; r_jj(kk) = 0; r_pos(kk) = 1.37e-10; r_neg(kk) = 2.58e-10 ;
%         dataloc = 'C:\Postdoctoral research\DataLit\Data_Single_Salts\';
%         datafile = strcat(dataloc,Solute,'.m');
%         run(datafile)

%         m_data = x;     % x = molality
%         o_data = Q;     % Q = Osmotic coefficient
%         aw_data = exp(-o_data.*v.*Mw.*m_data);
%         x_data = (v.*m_data)./(v.*m_data+1/Mw);
%         root_x_data = sqrt(x_data);

    elseif isequal(Solute,'Na2SO4')
        j(kk) = 15;  zpos(kk) = 1; zneg(kk) = 2; v(kk) = 3; rho(kk) = 8.110; n(kk) = 3; mu_j(kk) = 47.02; r_jw(kk) = 7.32e-10; r_jj(kk) = 4.3448e-10; r_pos(kk) = 1.02e-10; r_neg(kk) = 2.58e-10 ;

%         dataloc = 'C:\Postdoctoral research\DataLit\Data_Single_Salts\';
%         datafile = strcat(dataloc,Solute,'.m');
%         run(datafile)

%         m_data = x;     % x = molality
%         o_data = Q;     % Q = Osmotic coefficient
%         aw_data = exp(-o_data.*v.*Mw.*m_data);
%         x_data = (v.*m_data)./(v.*m_data+1/Mw);
%         root_x_data = sqrt(x_data);

    elseif isequal(Solute,'H2CO3')
        j(kk) = 16;  zpos(kk) = 1; zneg(kk) = 2; v(kk) = 3; rho(kk) = 13.00; n(kk) = 2; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 4.02E-10; r_pos(kk) = 0.21E-10; r_neg(kk) = 1.79E-10 ;
%     elseif isequal(Solute,'NH4HC2O4')
%         j(kk) = 13;  zpos(kk) = 1; zneg(kk) = 1; v(kk) = 2; rho(kk) = 13.00; n(kk) = 2; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 4.02E-10; r_pos(kk) = 1.78E-10; r_neg(kk) = 1.56E-10 ;
%     elseif isequal(Solute,'(NH4)2C2O4')
%         j(kk) = 14;  zpos(kk) = 1; zneg(kk) = 2; v(kk) = 3; rho(kk) = 13.00; n(kk) = 2; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 4.02E-10; r_pos(kk) = 1.78E-10; r_neg(kk) = 1.79E-10 ;
%     elseif isequal(Solute,'H2C2O4')
%         j(kk) = 15;  zpos(kk) = 1; zneg(kk) = 2; v(kk) = 3; rho(kk) = 13.00; n(kk) = 2; mu_j(kk) = 0; r_jw(kk) = 0; r_jj(kk) = 4.02E-10; r_pos(kk) = 0.21E-10; r_neg(kk) = 1.79E-10 ;
        
%         elseif isequal(Solute(kk),{'none'})
%         j(kk) = 1; zpos(kk) = 1; zneg(kk) = 1; v(kk) = 1; rho(kk) = 1; n(kk) = 2; mu_j(kk) = 1; r_jw(kk) = 19e-10; r_jj(kk) = 0; r_pos(kk) = 1e-10; r_neg(kk) = 1e-10 ;
    end
    
    if j(kk) >= 6 && j(kk) <= 13
        mu_j(kk) = q1*e*((r_pos(kk)+r_neg(kk))/2 + r_jj(kk))/D;
        r_jw(kk) = r_jj(kk) + (r_pos(kk)+r_neg(kk))/2 + r_ww/2;
% %     elseif j(kk) > 13
% %         mu_j(kk) = 0.866*q2*e*((r_pos(kk)+r_neg(kk)) + r_jj(kk))/D;
% %         r_jw(kk) = 0.866*(r_jj(kk) + (r_pos(kk)+r_neg(kk))) + r_ww/2;
    end
    
    for ii = 1:n(kk)-1
        deltaE(kk,ii) = ((mu_j(kk)*mu_w*D^2)/((1.113e-10)*(r_jw(kk)+(ii-1)*r_ww)^3))-((mu_w^2*D^2)/((1.113e-10)*(ii*r_ww)^3));
        C(kk,ii) = exp(deltaE(kk,ii)/(k*T)); % Energy parameter
        % C(kk,ii) = 1;
    end
     
end


for kk = 1:N-1
%     prompt = 'Enter ratio of current solute to last solute in the form "number" as in (number:1)';
%     soluteratio(kk) = input(prompt);

    soluteratio(kk,:) = molarRatio(kk);
end

assignin('base', 'soluteratio', soluteratio);

t = length(aw);

for i = 1:t
    m(N,i) = NR(j,rho,v,n,zpos,zneg,C,aw,N,soluteratio);
%     m1(N,i)= m(N,i);
    for kk = 1:N-1
        m(kk,i) = m(N,i).*soluteratio(kk,:);
%         m2(kk,i)= m(kk,i);
    end
    
    m_total(i) = 0;
    for k = 1:N
        m_total(i) = m_total(i) + m(k,i)*v(k);
    end
%     o(i) = (-log(aw(i)))./(Mw*m_total(i));
%     root_m(i) = power(m_total(i),0.5);
%     x (i) = m_total(i)./(m_total(i)+1/Mw); % Mole fraction of solute
%     root_x(i) = sqrt(x(i));
    m_dry_total(i) = sum(m(:,i));
    
%     for k = 1:N
%         x_each_solute (k,i) = (v(k).*m(k,i))./((m_total(i))+1/Mw);
% %         x_each_solute (k,i) = (m(k,i))./((v(k).*m(k,i))+1/Mw);
%     end
        
    for kkk = 1:N
        xstar(kkk) = m(kkk,i)./m_dry_total(i);
    end
    for p =1:N
        [a(p,i),Ix(p,i),f_Inf(p,i),gamma_Inf(p,i),K(p,i)] = SoluteActivity(N,i,j(p),rho,v,n(p),zpos,zneg,C(p,:),m,aw,p,xstar(p));
        if j(p) >= 6 && j(p) <= 15
            f_FS(p,i) = (m_total(i)+(1./Mw)).*((a(p,i)./K(p,i))./(((zpos(p).*m(p,i)).^zpos(p)).*((zneg(p).*m(p,i)).^zneg(p)))).^(1./v(p));  % Fused salt reference state; Mole fraction basis
        else
            f_FS(p,i) = (m_total(i)+(1./Mw)).*((a(p,i)./K(p,i))./m(p,i)); % Fused salt reference state; Mole fraction basis
        end
%         Inf_ratio(p,i) = gamma_Inf(p,i)./f_Inf(p,i);
    end
end

% aw=aw';
% assignin('base', 'WaterActivity', aw);
% m=m';
% 
% assignin('base', 'SoluteMolality', m);
% a=a';
% x=x';
% f_FS=f_FS';
% f_Inf=f_Inf';
% gamma_Inf=gamma_Inf';
% % % % K=K';
% % % % Ix=Ix';
% % % % IxFS=IxFS';
% % % % Kw=Kw';
% % % % awbar=awbar';
% % % % K_Inf=K_Inf';
% % % % root_x=root_x';
% % % % xstar=xstar';
% % % % x_each_solute=x_each_solute';
% % % % xpos=xpos';
% % % % xneg=xneg';

% Inf_ratio=Inf_ratio';

% mu_j;
% m1=m1';
% m2=m2';

% assignin('base', 'MixtureMolality', m_total);
% assignin('base', 'OsmoticCoeff', o);
% assignin('base', 'SoluteActivity', a);
% assignin('base', 'MoleFraction', x);
% assignin('base', 'ActCoeff_f_FS', f_FS);
% assignin('base', 'ActCoeff_f_Inf', f_Inf);
% assignin('base', 'gamma_Inf', gamma_Inf);
% assignin('base', 'xstar', xstar);

end

function [mout] = NR(j,rho,v,n,zpos,zneg,C,aw,N,soluteratio)

dm = 0.0000000001; % step size for derivative of function f
m1 = 1; % Give an initial starting molality guess.

% Newton raphson method

f = @(m) molalCalc(j,rho,v,n,zpos,zneg,C,m,aw,N,soluteratio) - m;
df = @(m)(f(m+dm)-f(m))/dm;
for i = 1:10
    m2 = abs(m1 - f(m1)/df(m1));
%     error(i) = (m1-m2).^2./m1.^2;
    m1 = m2;
end

mout = m1;

% assignin('base', 'ErrorMolality', error);

end

function [OUT] = molalCalc(j,rho,v,n,zpos,zneg,C,m,aw,N,soluteratio)

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

% disp(Ix); fprintf('\b'); disp(Kw); fprintf('\b'); disp(awbar);

% Calculate molality Using EQN.27
for p =1:N
    
%     if j(p) >= 1 && j(p) <= 5
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
%     else
%         NumorSum(p)=0;
%         for k=1:(n(p)-1)
%             NumorSum(p)=NumorSum(p)+(awbar.^k).*(1-C(p,k)).*prod(C(p,1:k-1));
%         end
%         DenomSum(p)=0;
%         for k=1:(n(p)-2)
%             DenomSum(p)=DenomSum(p)+k*(awbar.^(k-1)).*prod(C(p,1:k));
%         end
%         Denom(p)=(1-awbar).^2.*(DenomSum(p))+(n(p)-1-(n(p)-2).*awbar).*awbar.^(n(p)-2).*prod(C(p,1:n(p)-1));
%         m0(p)=((1-awbar)/(Mw*v(p)*awbar))*(1-NumorSum(p))/Denom(p);
%     end
end

assignin('base', 'BinaryMolality', m0);

mixturemodel=1./m0(N);
for k=1:N-1
    mixturemodel=mixturemodel+soluteratio(k)/m0(k);
end

OUT = 1./mixturemodel;

end

function [a,Ix,f_Inf,gamma_Inf,K] = SoluteActivity(N,i,j,rho,v,n,zpos,zneg,C,m,aw,p,xstar)

Mw=0.0180152; % Molecular weight water (kg/mol)
Ax=2.917; % Debye-Huckel coeff. mol frac @ 298.15K

IxTop=0;
for k=1:N
    IxTop=IxTop+m(k,i)*zpos(k)*zneg(k).*v(k);
end

IxBottom=0;
for k=1:N
    IxBottom=IxBottom+m(k,i).*v(k);
end

Ix = (1/2)* (IxTop)./(IxBottom+1/Mw);

% Ixref = (1/2)* zpos(p).*zneg(p);
% IxFS = (1/2)* zpos(p).*zneg(p);
Ixref = (1/2)* (IxTop)./(IxBottom);
IxFS = (1/2)* (IxTop)./(IxBottom);

KwTop=0;
for k=1:N
    KwTop=KwTop+(m(k,i)*zpos(k)*zneg(k).*v(k))./(1+rho(k)*Ix.^0.5);
end

Kw = exp((Ax.*(Ix.^0.5).*(KwTop))./(IxBottom+1/Mw));
awbar=aw./(Kw);

%
% IxTop = m(1,i)*zpos(1)*zneg(1).*v(1) + m(2,i)*zpos(2)*zneg(2).*v(2);
% IxBottom = m(1,i)*v(1) + m(2,i)*v(2);
% KwTop = (m(1,i)*zpos(1)*zneg(1).*v(1))./(1+rho(1)*Ix.^0.5) + (m(2,i)*zpos(2)*zneg(2)*v(2))./(1+rho(2)*Ix.^0.5);

% disp(Ix); fprintf('\b'); disp(Kw); fprintf('\b'); disp(awbar);

DenSum=0;
for k=1:(n-1)
    DenSum=DenSum+(awbar.^k).*(1-C(k)).*prod(C(1:k-1));
end

abar0 = ((1-awbar)./(1-DenSum)).^v(p);

if j >= 1 && j <= 5
    
%     K = exp((Ax*Ix^0.5*KwTop)/(IxBottom+1/Mw)); % mixture with only neutral species
    K=exp((-Ax.*(IxFS-Ix).*KwTop)./((IxBottom+1./Mw).*Ix.^0.5)); % mixture with atleast one eletrolyte

else
    SumKjTerm=0;
    for k=1:N
        SumKjTerm=SumKjTerm+(m(k,i)*zpos(k)*zneg(k).*v(k))./(2.*Ix.^0.5.*(1+rho(k)*Ix.^0.5));    % original
%         SumKjTerm=SumKjTerm+(m(k,i).*v(k).*Ixref)./(Ix.^0.5.*(1+rho(k)*Ix.^0.5));   % trying IxFS and Ixref in place of zz/2.  
    end
    
%     K = (exp(-zpos(p).*zneg(p)*Ax*(2/rho(p)*log((1+rho(p)*Ix^0.5)/(1+rho(p)*Ixref^0.5))+((1-((2*Ix)./(zpos(p).*zneg(p))))/(IxBottom+1/Mw))*SumKjTerm)))^v(p);   % Eq. 16 from Paper 3
    K = (exp(-zpos(p).*zneg(p)*Ax*(2/rho(p)*log((1+rho(p)*Ix^0.5)/(1+rho(p)*Ixref^0.5))+((1-Ix./IxFS)/(IxBottom+1/Mw))*SumKjTerm)))^v(p);
end


a = xstar.*K.*abar0;

% Reference state --> Infinite dilution
Ixref = 0;

for k=1:(n-1)
    Cprod=prod(C(1:k));
end

if j >= 1 && j <= 5
    
%     K_Inf = exp((Ax*Ix^0.5*KwTop)/(IxBottom+1/Mw)); % mixture with only neutral species
    K_Inf=exp((-Ax.*(IxFS-Ix).*KwTop)./((IxBottom+1./Mw).*Ix.^0.5)); % mixture with atleast one eletrolyte

    f_Inf = ((IxBottom+(1./Mw))./IxBottom).*abar0.*Cprod.*K_Inf; % Infinite dilution reference state; Mole fraction basis Activity Coefficient for Non-dissociating species
else
%     K_Inf = (exp(-zpos(p).*zneg(p)*Ax*(2/rho(p)*log((1+rho(p)*Ix^0.5)/(1+rho(p)*Ixref^0.5))+((1-((2*Ix)./(zpos(p).*zneg(p))))/(IxBottom+1/Mw))*SumKjTerm)))^v(p);   % Eq. 16 from Paper 3
    K_Inf = (exp(-zpos(p).*zneg(p)*Ax*(2/rho(p)*log((1+rho(p)*Ix^0.5)/(1+rho(p)*Ixref^0.5))+((1-Ix./IxFS)/(IxBottom+1/Mw))*SumKjTerm)))^v(p);
    f_Inf = ((IxBottom+(1./Mw))./IxBottom).*((abar0.*K_Inf).^(1./v(p))).*Cprod;  % Infinite dilution reference state; Mole fraction basis Activity Coefficient for Dissociating species
end

% Molality basis Inf dilution reference state Activity coefficient (eq. 30)
gamma_Inf = (Cprod./(Mw.*IxBottom)).*((abar0.*K_Inf).^(1./v(p)));


% Reference state --> Fused salt

xpos = (zneg(p).*m(p,i))./(IxBottom+(1/Mw));
xneg = (zpos(p).*m(p,i))./(IxBottom+(1/Mw));

end
