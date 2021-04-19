function [result,s1,s2] = Mie_tetascan(m, x, nsteps)

% x=12.9559; m=1.37486335230751 + 0.00500000000000000i; nsteps=19;


% Computation and plot of Mie Power Scattering function for given 

% complex refractive-index ratio m=m'+im", size parameters x=k0*a, 

% according to Bohren and Huffman (1983) BEWI:TDD122

% C. M�tzler, May 2002.



nsteps=nsteps;

m1=real(m); m2=imag(m);

nx=(1:nsteps); dteta=pi/(nsteps-1);

teta=(nx-1).*dteta;

    for j = 1:nsteps

        u=cos(teta(j));

        a(:,j)=Mie_S12(m,x,u);
        s1(j) = a(1,j);
        s2(j) = a(2,j);

        SL(j)= real(a(1,j)'*a(1,j));

        SR(j)= real(a(2,j)'*a(2,j));

    end

y=[teta teta+pi;SL SR(nsteps:-1:1)]';



% % % % polar(y(:,1),y(:,2))
% % % % 
% % % % title(sprintf('Mie angular scattering: m=%g+%gi, x=%g',m1,m2,x));
% % % % 
% % % % xlabel('Scattering Angle')

result=y; 

