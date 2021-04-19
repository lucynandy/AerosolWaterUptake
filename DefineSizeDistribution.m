% Random Parameter Generation by Mehdi Rezaeianzadeh
% Uniform Distribution  
% Log Normal Distribution
% Triangular Distribution

% % % % clc
% % % % clear all

function [r_lognorm] = DefineSizeDistribution(n,m,k,mu,sigma)

%fid = fopen('input1.txt');
%[m,n,Xl,Xu,Xll,Xuu,Rl,Ru,a,b,c] = textread('input1.txt','%f %f %f %f %f %f %f %f %f %f %f')
%mm=input('Number of Parameters:');
% % % % [n] = textread('input_control_sims.txt','n=%f %*[^\n]',1,'whitespace','\n'); %number of parameters and is fixed here
                                      
% [m,k] = textread('input_control.txt','m=%f k=%f %*[^\n]',4,'whitespace','\n'); %The users need to set hom many parameters they have
% % % % [m,k] = textread('input_control.txt','m=%f k=%f %*[^\n]',1,'whitespace','\n'); %The users need to set hom many parameters they have
                                   
% input_control=xlsread('input_control');   
% m=input_control(:,1);
% n=input_control(:,2);
% k=input_control(:,3);

%m=input('m= Enter a number from 1 to 3 to choose among "Uniform=1", "Log Normal=2", "Triangular=3" Distributions:')
%n=input('n= The number of random variables:')
for f=1:length(m(:,1))
switch m(f)
    case 1  %Random data generation based on "Uniform Distribution"
%         Xl=input('Min of the parameter range for Uniform Distribution:')
%         Xu=input('Max of the parameter range for Uniform Distribution:')
f1= sum(m == 1); % counts the number of m=1's
         [Xl,Xu] = textread('input_uniform_dist.txt','Xl=%f Xu=%f %*[^\n]',f1,'whitespace','\n');
                                      
% input_uniform_dist=xlsread('input_uniform_dist');
% Xl=input_uniform_dist(:,1);
% Xu=input_uniform_dist(:,2);

f1= sum(m == 1);
r_uniform=zeros(n(1),f1);
        for ii=1:f1
        r_uniform(:,ii)= (Xu(ii)-Xl(ii)).*rand(n(1),1) + Xl(ii);
        end
%         output= fopen('output_uniform.txt','w');
% %         fprintf(output,'%6s %6s %6s %6s\r\n','r_uniform');
%         fprintf(output,'%6f \r\n',r_uniform);
        dlmwrite('output_uniform.txt',r_uniform,'delimiter','\t','precision',6)
      %xlswrite('r_uniform.xlsx',r_uniform)
    case 2  %Random data generation based on "Log Normal Distribution"
%         Xl=input('Min of the parameter range for Log Normal Distribution:')
%         Xu=input('Max of the parameter range for Log Normal Distribution:')
%         Rl=input('Lower percentile for Log Normal Distribution:')
%         Ru=input('Upper percentile for Log Normal Distribution:')

          switch k(f) % to select between two scenarios of having "mu" and "sigma" or not
              case 1 
                  f21= sum(m == 2 & k==1);
% % % %                    [mu,sigma] = textread('input_lognorm_dist_mu.txt','mu=%f sigma=%f %*[^\n]',f21,'whitespace','\n');
% input_lognorm_dist_mu=xlsread('input_lognorm_dist_mu');
% mu=input_lognorm_dist_mu(:,1);
% sigma=input_lognorm_dist_mu(:,2);

                  for iii=1:f21
                  r_lognorm(:,iii)= logninv(rand(n(1),1),mu(iii),sigma(iii));  %Lognormal inverse cumulative distribution function
                  end
%                   output = fopen('output_logN_mu.txt','w');
%                   fprintf(output,'%6s \r\n','r_lognorm_mu');
%                   fprintf(output,'%6f \r\n',r_lognorm_mu);
% % % %                   dlmwrite('output_lognorm_mu.txt',r_lognorm_mu,'delimiter','\t','precision',6)
              case 2
                  f22= sum(m == 2 & k==2);
                  Rl=0.1;    Ru=99.9;
                  Xl=mu;    Xu=sigma;
                  % % % %                    [Xl,Xu,Rl,Ru] = textread('input_lognorm_dist.txt','Xl=%f Xu=%f Rl=%f Ru=%f %*[^\n]',f22, ...
% % % %                                        'whitespace','\n');
% input_lognorm_dist=xlsread('input_lognorm_dist');
% Xl=input_lognorm_dist(:,1);
% Xu=input_lognorm_dist(:,2);
% Rl=input_lognorm_dist(:,3);
% Ru=input_lognorm_dist(:,4);

for iiii=1:f22
                  Finv_Rl(iiii)=norminv(Rl(iiii)/100,0,1);%Normal inverse cumulative distribution function
                  Finv_Ru(iiii)=norminv(Ru(iiii)/100,0,1);
                  mu(iiii)=(log(Xu(iiii))*Finv_Rl(iiii)-log(Xl(iiii))*Finv_Ru(iiii))/(Finv_Rl(iiii)-Finv_Ru(iiii));% Mean calculation
                  sigma(iiii)=log(Xu(iiii)/Xl(iiii))/(Finv_Ru(iiii)-Finv_Rl(iiii));% StD calculation
                
                  r_lognorm(:,iiii)= logninv(rand(n(1),1),mu(iiii),sigma(iiii));
end
                  %Lognormal inverse cumulative distribution function
%                   output = fopen('output_LogN.txt','w');
%                   fprintf(output,'%6s \r\n','r_lognorm');
%                   fprintf(output,'%6f \r\n',r_lognorm);
% % % %                   dlmwrite('output_lognorm.txt',r_lognorm,'delimiter','\t','precision',6)
          end
    case 3 %Random data generation based on "Triangular Distribution"
%         a=input('input "a" as lower limit:');%lower limit
%         c=input('input "c" as Peak Location:');%Peak Location
%         b=input('input "b" as Upper limit:');%Upper limit
         f3= sum(m == 3);
         [a,c,b] = textread('input_triangular_dist.txt','a=%f c=%f b=%f %*[^\n]',f3, ...
                                       'whitespace','\n');
                                       
%   input_triangular_dist=xlsread('input_triangular_dist');
%   a=input_triangular_dist(:,1);
%   c=input_triangular_dist(:,2);
%   b=input_triangular_dist(:,3);
 
        U=rand(n(1),1);
        Fc(:,1)=(c(1:f3,1)-a(1:f3,1))./(b(1:f3,1)-a(1:f3,1));
 
        for i=1:length(U)
            if U(i)<Fc(1:f3)
               r_triangular(i,1:f3)=a(1:f3)+sqrt(U(i).*(b(1:f3)-a(1:f3)).*(c(1:f3)-a(1:f3)));
            else
               r_triangular(i,1:f3)=b(1:f3)-sqrt((1-U(i)).*(b(1:f3)-a(1:f3)).*(b(1:f3)-c(1:f3)));
            end
        end
%         output = fopen('output_triangular.txt','w');
%         fprintf(output,'%6s \r\n','r_triangular');
%         fprintf(output,'%6f \r\n',r_triangular);
             dlmwrite('output_triangular.txt',r_triangular,'delimiter','\t','precision',6)
end
end
