% Convergence approximation discontinuous functions
% Modified 
% piecewise analytic function

clear all


fun = @(x) (0.005.*x.^4.-2.*x).*(-5 <= x & x < -2) + (atan(exp(x.^2))).*(-2 <= x & x < 1) + (log(sin(x))).*(1 <= x & x < 3) + (acosh(x.^3)).*(3 <= x & x <= 5);
tol = 1e-13; % tolerance 1e-13/1e-14


M = 5000; % Sample points per domain
Z = linspace(-5,5,M);


G = linspace(-5,5,5000)';

d = [abs((atan(exp((-2).^2)))-(0.005.*(-2).^4.-2.*(-2)));
     abs(log(sin(1)) - (atan(exp(1.^2))));
     abs((acosh(3.^3)) - (log(sin(3))))];   % Jumps

S = @(x) sGibbs([-2 1 3],d,10,x); % map 

AllmaxErrorAAA = [];


for m = 1:100%max type (m-1,m-1)
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),S(Z),tol,m); 
rS = @(x) r(S(x));

errfunAAA = @(x) abs(fun(x)-rS(x));

AllmaxErrorAAA = [AllmaxErrorAAA, max(abs(fun(G)-rS(G)))];
end



figure 
plot(AllmaxErrorAAA,'o-','Color',[0 0 1],'LineWidth',2)
hold on
leg{1} = 'AAA Approximant of type(m-1,m-1)';
title_string = strcat(['Appr. discontinuous functions, Modified, piecewise analytic function']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','LaTex')
ylabel('max error','Interpreter','LaTex')
set(gca, 'YScale', 'log')
grid on
hold off

% Real part of AAA Approximant
rreal = @(x) real(rS(x)); 

% plot with disconnectet Z without Gibbs pheno.
figure 
fplot(fun,[-5,5],'-','Color',[1 0 0],'LineWidth',3) % Function
hold on 
title_string = strcat(['Appr. discontinuous functions, Modified, piecewise analytic function']);
title(title_string,'Interpreter','LaTex','FontSize',20);
fplot(rreal,[-5,5],'-','Color',[0 0 1],'LineWidth',2) % AAA-App
leg{1} = 'Function';
leg{2} = ['AAA Approximant of type (' num2str(length(pol)) ',' num2str(length(pol)) ')'];
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('Re$(z)$','Interpreter','LaTex')
ylabel('Re$(r(z))$','Interpreter','LaTex')
grid on
hold off