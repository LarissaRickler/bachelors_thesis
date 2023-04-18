% Convergence of application 5
% 5 Appr. in connected domains
% 1/bessel

clear all

M = 2000; % # sample points, paper 2000
fun = @(x) 1./besselj(0,x); %1/J_0 with J_0 bessel fct
tol = 1e-14; % 1e-14;  % 1e-13/1e-14

% Sample points
% M random points in complex plane (0,10)+i(-1,1)
Z = 10*rand(M,1) + 1i*(-1+2*rand(M,1));                    

% Grid 200x1000 points in complex plane (0,10)+i(-1,1)
[GX,GY] = meshgrid(linspace(0,10,1000),linspace(-1,1,200));
G1 = GX + 1i*GY;
G = G1(:)';




AllmaxErrorAAA = [];

for m=1:40 %max type (m-1,m-1), in paper m = 15
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,tol,m); 

errfunAAA = @(x) abs(fun(x)-r(x));

AllmaxErrorAAA = [AllmaxErrorAAA, max(errfunAAA(G))];
end


figure 
plot(AllmaxErrorAAA,'o-','Color',[0 0 1],'LineWidth',2)
hold on
leg{1} = 'AAA Approximant of type(m-1,m-1)';
title_string = strcat('Appr. in connected domains, $f(z) = \frac{1}{J_0}$ with  bessel func. $J_0$');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','LaTex')
ylabel('max error','Interpreter','LaTex')
set(gca, 'YScale', 'log')
grid on
hold off

