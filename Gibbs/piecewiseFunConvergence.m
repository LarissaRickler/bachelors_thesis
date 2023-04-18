% Convergence approximation discontinuous functions
% piecewise analytic function

clear all


fun = @(x) (0.005.*x.^4.-2.*x).*(-5 <= x & x < -2) + (atan(exp(x.^2))).*(-2 <= x & x < 1) + (log(sin(x))).*(1 <= x & x < 3) + (acosh(x.^3)).*(3 <= x & x <= 5);
tol = 1e-13; % tolerance 1e-13/1e-14


M = 5000; % Sample points per domain
Z = linspace(-5,5,M);


G = linspace(-5,5,5000);




AllmaxErrorAAA = [];


for m = 1:100 %max type (m-1,m-1)
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,tol,m); 

errfunAAA = @(x) abs(fun(x)-r(x));

AllmaxErrorAAA = [AllmaxErrorAAA, max(errfunAAA(G))];
end


figure 
plot(AllmaxErrorAAA,'o-','Color',[0 0 1],'LineWidth',2)
hold on
leg{1} = 'AAA Approximant of type(m-1,m-1)';
title_string = strcat(['Appr. discontinuous functions, Sawtooth']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','LaTex')
ylabel('max error','Interpreter','LaTex')
set(gca, 'YScale', 'log')
grid on
hold off

% Real part of AAA Approximant
rreal = @(x) real(r(x)); 

% plot with disconnectet Z without Gibbs pheno.
figure 
fplot(fun,[-5,5],'-','Color',[1 0 0],'LineWidth',3) % Function
hold on 
title_string = strcat(['Appr. discontinuous functions, Sawtooth']);
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