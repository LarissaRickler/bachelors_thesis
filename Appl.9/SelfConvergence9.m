% Convergence of application 9
% 9 ISS model form Chahlaoui and Van Dooren
% fun(z) = c'(zI-A)^(-1)b (a rat. fun. of type(348,348) 



clear all

load('ISS.mat') % loads A,B,C of model


dim = 270;
tol = 1e-5; % tolerance 1e-13/1e-14


% Sample points
M = 500; % % 2*M sample points, paper n = 500
Z1 = 1i*logspace(-2,2,M); % M logarithmically spaced points from 10?2 i to 102 i
Z2 = -Z1; % complex conjugates
Z = [Z1 Z2]';     

% Grid
G = 1i*linspace(-100,100,20000)';


% load funZ: evaluation of fun at Z
load('funZ9.mat')

% load funG: evaluation of fun at G
load('funG9.mat')

it = 1;

figure

for in = 1:3 % row of C
    for out = 1:3 % column of B
funZ = funZ9(3*(in-1)+out,:)';
funG = funG9(3*(in-1)+out,:)';


AllmaxErrorAAA = [];
for m = 1:120 % type (m-1,m-1), paper between 0 and 46
[r,pol,res,zer,z,f,w,errvec] = aaa(funZ,Z,tol,m); 

rG = r(G);
err = abs(funG-rG);

AllmaxErrorAAA = [AllmaxErrorAAA, max(err)];
end


subplot(3,3,it)
plot(AllmaxErrorAAA,'o-','Color',[0 0 1],'LineWidth',2)
hold on
leg{1} = 'AAA Approximant of type (m-1,m-1)';
title_string = strcat(['$i = $' num2str(in) ', $j = $' num2str(out)]);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','LaTex')
ylabel('max error','Interpreter','LaTex')
set(gca, 'YScale', 'log')
grid on
hold off



% subplot(3,3,it)
% plot(imag(G),abs(funG),'-','Color',[1 0 0],'LineWidth',3)
% hold on 
% plot(imag(G),abs(rG),'-','Color',[0 0 1],'LineWidth',2) % AAA-App
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% xlim([-100 100])
% leg{1} = 'Function';
% leg{2} = ['AAA Approximant of type (' num2str(m-1) ', ' num2str(m-1) ')'];
% title_string = strcat(['ISS model: input ' num2str(in) ', output ' num2str(out)]);
% title(title_string,'Interpreter','LaTex','FontSize',14);
% legend(leg,'Interpreter','LaTex','FontSize',14,'Location','EastOutside');
% set(gca,'FontSize',14);
% xlabel('Im($z$)','Interpreter','LaTex')
% ylabel('$|f(z)|$ or $|r(z)|$','Interpreter','LaTex')
% grid on 
% hold off



it = it+1;
    end 
end
