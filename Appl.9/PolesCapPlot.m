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

it = 2;
figure

for in = 1 % row of C
    for out = 1 % column of B
funZ = funZ9(3*(in-1)+out,:)';
funG = funG9(3*(in-1)+out,:)';


subplot(3,2,1)
plot(imag(G),abs(funG),'-','Color',[1 0 0],'LineWidth',1)
hold on 
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
title_string = strcat(['function $f_{1,1}$']);
title(title_string,'Interpreter','LaTex','FontSize',14);
set(gca,'FontSize',14);
xlabel('Im($z$)','Interpreter','LaTex')
ylabel('$|f(z)|$','Interpreter','LaTex')
ylim([1e-7 1])
xlim([0.01 100])
grid on 
hold off


for m = 3:2:12 % type (m-1,m-1), paper between 0 and 46
[r,pol,res,zer,z,f,w,errvec] = aaa(funZ,Z,tol,m); 
rG = r(G);

subplot(3,2,it)
plot(imag(G),abs(rG),'-','Color',[0 0 1],'LineWidth',1) % AAA-App
hold on 
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
title_string = strcat(['approximant $r$ of type $($' num2str(m-1) '$,$' num2str(m-1) ')']);
title(title_string,'Interpreter','LaTex','FontSize',14);
set(gca,'FontSize',14);
xlabel('Im($z$)','Interpreter','LaTex')
ylabel('$|r(z)|$','Interpreter','LaTex')
ylim([1e-7 1])
xlim([0.01 100])
grid on 
hold off

it = it+1;
end



    end 
end
