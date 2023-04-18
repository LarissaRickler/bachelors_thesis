
fun = @(x) 1./(1+25.*x.^2);
M = 10; % # sample points, paper 1000
tol = 1e-14; % 1e-13/1e-14

Z = linspace(-1,1,M);
ZR = linspace(-1,1,M/2);

mr = 6; % type(mr-1,mr-1)
mp = 2*(mr-1); % degree mp

[r,pol,res,zer,z,f,w,errvec] = aaa(fun(ZR),ZR,tol,mr);
% [r,arat,brat,murat,nurat,polesrat,residuesrat] = ratdisk(fun,mr-1,mr-1,M,tol);
rreal = @(x) real(r(x));

pc = polyfit(Z,fun(Z),mp);
p = @(x) polyval(pc,x);


figure 
fplot(fun,[-1,1],'-','Color',[1 0 0],'LineWidth',5) % fun
hold on 
fplot(rreal,[-1,1],'--','Color',[0 0 1],'LineWidth',3) % AAA-app
fplot(p,[-1,1],'--','Color',[1 0 1],'LineWidth',3) % Poly


leg{1} = 'Runge function';
leg{2} = sprintf('Rational interpolant of type (%i,%i)', mr-1 ,mr-1);
leg{3} = sprintf('Polynomial interpolant of degree %i', mp);
title_string = strcat('Runge function');
% title(title_string,'Interpreter','LaTex','FontSize',12);
legend(leg,'Interpreter','LaTex','FontSize',14,'Location','EastOutside');
set(gca,'FontSize',14);
set(0,'defaulttextInterpreter','latex')
grid on
hold off
