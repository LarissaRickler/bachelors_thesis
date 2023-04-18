%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   Gibbs phenomenon  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % % % % % Example paper % % % % %

clear all

% function
g = @(x) (exp(1./(x + 5.5)).*(-5 <= x & x < -3)) + (cos(3.*x)).*(-3 <= x & x < 2) + (-x.^3./30 + 2).*(2 <= x & x <= 5);


% interval
a = -5; b = 5;

m = 100;

% AAA-interpolant
nAAA = 10000; % number of sample points
xAAA = linspace(a,b,nAAA); % equispaced sample points
F = g(xAAA); % function values at sample points
[r,pol,res,zer,z,f,w,errvec] = aaa(g,xAAA,1e-13,m); % Interpolant 

% AAA-Interpolant with fake nodes
nAAAfn = 10000; % number of sample points
xAAAfn = linspace(a,b,nAAAfn); % equispaced sample points
k = 50; % shifting parameter
S = @(x) sGibbs([-3, 2],[2.40295496, 0.7731630467],k,x); % fake nodes map
Gfn = g(xAAAfn); % function values at sample points
[rfn,polfn,resfn,zerfn,zfn,ffn,wfn,errvecfn] = aaa_mod(Gfn,xAAAfn,S,1e-13,m); % interpolant

% plot
figure
fplot(g, [a, b], '-', 'Color', [1,0,0], 'LineWidth', 3)
hold on
fplot(r, [a, b], '--', 'Color', [0,0,1], 'LineWidth', 3)
fplot(rfn, [a, b], '--', 'Color', [1,0,1], 'LineWidth', 3)

leg{1} = 'function';
leg{2} = 'rational Interpolant';
leg{3} = 'rational Interpolant with fake nodes';
title_string = strcat('Sample');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off
