% Convergence approximation discontinuous functions
% Modified
% Square wave

clear all


fun = @(x) sign(sin((2.*pi.*real(x))./2));
tol = 1e-13; % tolerance 1e-13/1e-14


M = 5000; % Sample points per domain
Z = linspace(0,5,M);
Z = Z(2:end);

G = linspace(0,5,5000);
G = setdiff(G,[0,1,2,3,4])';

S = @(x) sGibbs([1 2 3 4],[2 2 2 2],10,x);

AllmaxErrorAAA = [];


for m = 1:100%max type (m-1,m-1)
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),S(Z),tol,m); 
rS = @(x) r(S(x));

errfunAAA = @(x) abs(fun(x)-rS(x));

AllmaxErrorAAA = [AllmaxErrorAAA, max(errfunAAA(G))];
end



figure 
plot(AllmaxErrorAAA,'o-','Color',[0 0 1],'LineWidth',2)
hold on
leg{1} = 'AAA Approximant of type(m-1,m-1)';
title_string = strcat(['Appr. discontinuous functions, Modified, $f(z) = $sign$(\sin(\frac{2\pi}{2}$Re$(z)))$']);
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
fplot(fun,[0,5],'-','Color',[1 0 0],'LineWidth',3) % Function
hold on 
title_string = strcat(['Appr. discontinuous functions, Modified, Square Wave']);
title(title_string,'Interpreter','LaTex','FontSize',20);
fplot(rreal,[0,5],'-','Color',[0 0 1],'LineWidth',2) % AAA-App
leg{1} = 'Function';
leg{2} = ['AAA Approximant of type (' num2str(length(pol)) ',' num2str(length(pol)) ')'];
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('Re$(z)$','Interpreter','LaTex')
ylabel('Re$(r(z))$','Interpreter','LaTex')
grid on
hold off