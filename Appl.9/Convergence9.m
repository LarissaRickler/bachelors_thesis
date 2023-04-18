% Convergence of application 9
% 9 Clamped beam model form Chahlaoui and Van Dooren
% fun(z) = c'(zI-A)^(-1)b (a rat. fun. of type(348,348) 



clear all

load('ClampedBeam.mat')


dim = 348;
fun = @(x) cAb(C,A,B,x);
tol = 1e-5; % tolerance 1e-13/1e-14



% Sample points
M = 500; % % 2*M sample points, paper n = 500
Z1 = 1i*logspace(-2,2,M); % M logarithmically spaced points from 10?2 i to 102 i
Z2 = -Z1; % complex conjugates
Z = [Z1 Z2]';     
% evaluation of fun at Z
funZ = fun(Z);

% Grid
G = 1i*linspace(-100,100,20000)';
funG = fun(G); % evaluation of fun at G




AllmaxErrorAAA = [];
for m = 1:50 % type (m-1,m-1), paper between 0 and 46
[r,pol,res,zer,z,f,w,errvec] = aaa(funZ,Z,tol,m); 

rG = r(G);
err = abs(funG-rG);

AllmaxErrorAAA = [AllmaxErrorAAA, max(err)];
end


figure 
plot(AllmaxErrorAAA,'o-','Color',[0 0 1],'LineWidth',2)
hold on
leg{1} = 'AAA Approximant of type (m-1,m-1)';
title_string = strcat(['Clamped beam model form Chahlaoui and Van Dooren']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','LaTex')
ylabel('max error','Interpreter','LaTex')
set(gca, 'YScale', 'log')
grid on
hold off


figure 
plot(imag(G),abs(funG),'-','Color',[1 0 0],'LineWidth',5)
hold on 
plot(imag(G),abs(rG),'-','Color',[0 0 1],'LineWidth',3) % AAA-App
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
leg{1} = 'Function';
leg{2} = ['AAA Approximant of type (' num2str(m-1) ', ' num2str(m-1) ')']';
title_string = strcat('Clamped beam model form Chahlaoui and Van Dooren');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off

