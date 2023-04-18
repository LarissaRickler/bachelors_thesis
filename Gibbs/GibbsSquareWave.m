% Gibbs square wave
% Figures


clear all
figure


fun = @(x) sign(sin((2.*pi.*real(x))./2)); % square wave
tol = 1e-13; % tolerance 1e-13/1e-14




M = 1000; % Sample points per circle

m = 80; % type (m-1,m-1)

Z = [];
for k = 1:5 % # domains
    Z = [Z; equisphere(k-0.5,0,0.2,M)];
end
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,tol,m); 

% Real part of AAA Approximant
rreal = @(x) real(r(x)); 

% plot with disconnectet Z without Gibbs pheno.
subplot(2,1,1)
fplot(fun,[0,5],'-','Color',[1 0 0],'LineWidth',3) % Function
hold on 
title_string = strcat(['Square Wave in Disconnected Domain']);
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





MG = 5000; % Sample points per circle

mG = 71; % type (m-1,m-1)

ZG = linspace(0,5,MG);
ZG = ZG(2:end);

[rG,polG,resG,zerG,zG,fG,wG,errvecG] = aaa(fun(ZG),ZG,tol,mG); 


% Real part of AAA Approximant
rrealG = @(x) real(rG(x)); 

% plot with disconnectet Z without Gibbs pheno.
subplot(2,1,2) 
fplot(fun,[0,5],'-','Color',[1 0 0],'LineWidth',3) % Function
hold on 
title_string = strcat(['Square Wave in Connected Domain']);
title(title_string,'Interpreter','LaTex','FontSize',20);
fplot(rrealG,[0,5],'-','Color',[0 0 1],'LineWidth',2) % AAA-App
leg{1} = 'Function';
leg{2} = ['AAA Approximant of type (' num2str(length(polG)) ',' num2str(length(polG)) ')'];
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('Re$(z)$','Interpreter','LaTex')
ylabel('Re$(r(z))$','Interpreter','LaTex')
grid on
hold off

