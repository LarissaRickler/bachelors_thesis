% Increasing depending on type
% Condition of problem and basis of application 5
% 5 Appr. in connected domains
% 1/bessel

clear all

M = 2000; % # sample points, paper 2000
fun = @(x) 1./besselj(0,x); %1/J_0 with J_0 bessel fct
tol = 1e-14; % 1e-14;  % 1e-13/1e-14

% Sample points
% M random points in complex plane (0,10)+i(-1,1)
Z = 10*rand(M,1) + 1i*(-1+2*rand(M,1));                    % "

% Grid 200x1000 points in complex plane (0,10)+i(-1,1)
[GX,GY] = meshgrid(linspace(0,10,1000),linspace(-1,1,200));
G1 = GX + 1i*GY;
G = G1(:)';



LambdaAll = [];
LambdaCAll = [];
LambdaCAllVF = [];
for m=1:70 %max type (m-1,m-1), in paper m = 15

[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,tol,m);

numL = @(x) sum(abs(w)./abs(x-z)); % numerator of Lebesgue fct
denomL = @(x) abs(sum(w./(x-z))); % denominator of Lebesgue fct
Lambdaf = numL(G)./denomL(G); % Lebesgue fct on points G
Lambda = max(Lambdaf); % "Lebesgue constant" 
LambdaAll = [LambdaAll, Lambda];

C = [];
Zm = setdiff(Z,z);
for k = 1:length(z) % Cauchy matrix AAA
C = [C 1./(Zm-z(k))];
end
LambdaC = cond(C, 2);
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
title_string = strcat(['Appr. in connected domains, $f(z) = \frac{1}{J_0}$ with  bessel func. $J_0$ ']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','Latex')
ylabel('Condition number $\Lambda^{(m)}$','Interpreter','Latex')
grid on
hold off

figure 
plot(LambdaCAll,'o-','Color',[0 0 1]) % AAA-App
hold on
leg{1} = 'Cauchy matrix AAA';
title_string = strcat(['Appr. in connected domains, $f(z) = \frac{1}{J_0}$ with  bessel func. $J_0$ ']);
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
title_string = strcat(['Appr. in connected domains, $f(z) = \frac{1}{J_0}$ with  bessel func. $J_0$ ']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
set(gca, 'YScale', 'log')
xlabel('$m$','Interpreter','Latex')
ylabel('Condition number $\kappa_{VF}$','Interpreter','Latex')
grid on
hold off