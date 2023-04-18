% Increasing depending on type
% Condition of problem and basis of application 9
% 9 Clamped beam model form Chahlaoui and Van Dooren
% fun(z) = c'(zI-A)^(-1)b (a rat. fun. of type(348,348)

clear all

load('ClampedBeam.mat')


dim = 348;
fun = @(x) cAb(C,A,B,x);
tol = 1e-14; % tolerance 1e-13/1e-14


% Sample points
M = 500; % % 2*M sample points, paper n = 500
Z1 = 1i*logspace(-2,2,M); % M logarithmically spaced points from 10?2 i to 102 i
Z2 = -Z1; % complex conjugates
Z = [Z1 Z2]';     
% evaluation of fun at Z
funZ = fun(Z);

% Grid
G = 1i*linspace(-100,100,20000) ;
                



LambdaAll = [];
LambdaCAll = [];
LambdaCAllVF = [];
for m = 1:50 %max type (m-1,m-1), in paper m = 15

[r,pol,res,zer,z,f,w,errvec] = aaa(funZ,Z,tol,m);

numL = @(x) sum(abs(w)./abs(x-z)); % numerator of Lebesgue fct
denomL = @(x) abs(sum(w./(x-z))); % denominator of Lebesgue fct
Lambdaf = numL(G)./denomL(G); % Lebesgue fct on points G
Lambda = max(Lambdaf); % "Lebesgue constant" 
LambdaAll = [LambdaAll, Lambda];

CAAA = [];
Zm = setdiff(Z,z);
for k = 1:length(z) % Cauchy matrix AAA
CAAA = [CAAA 1./(Zm-z(k))];
end
LambdaC = cond(CAAA, 2);
LambdaCAll = [LambdaCAll, LambdaC];


% 'Vector fitting'
CVF = [];
ZmVF = setdiff(Z,pol);
for k = 1:length(pol) % Cauchy matrix VF
CVF = [CVF 1./(ZmVF-pol(k))];
end
LambdaCVF = cond(CVF, 2);
LambdaCAllVF = [LambdaCAllVF, LambdaCVF];
end


figure 
plot(LambdaAll,'o-','Color',[1 0 0]) % Lambda
hold on
leg{1} = 'Condition of Problem';
title_string = strcat(['Clamped beam model form Chahlaoui and Van Dooren']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
set(gca, 'YScale', 'log')
xlabel('$m$','Interpreter','Latex')
ylabel('Condition number $\Lambda^{(m)}$','Interpreter','Latex')
grid on
hold off

figure 
plot(LambdaCAll,'o-','Color',[0 0 1]) % AAA-App
hold on
leg{1} = 'Cauchy matrix AAA';
title_string = strcat(['Clamped beam model form Chahlaoui and Van Dooren']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','Latex')
ylabel('Condition number $\kappa_{AAA}$','Interpreter','Latex')
grid on
hold off

figure
plot(LambdaCAllVF,'o-','Color',[0 1 1]) % VF
hold on
leg{1} = 'Cauchy matrix VF';
title_string = strcat(['Clamped beam model form Chahlaoui and Van Dooren']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
set(gca, 'YScale', 'log')
xlabel('$m$','Interpreter','Latex')
ylabel('Condition number $\kappa_{VF}$','Interpreter','Latex')
grid on
hold off