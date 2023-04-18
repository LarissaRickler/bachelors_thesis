% Increasing depending on type
% Condition of problem and basis of application 2 
% 2 Analytic fct with nearby brance points
% log(1.1 - x)

clear all


fun = @(x) log(1.1 -x);
M = 257; % # sample points, paper 256
tol = 1e-14; % 1e-13/1e-14

% Grid for Lambda
G1 = disk_grid(100,1,[0,0],1); % equispaced points in complex unit disk
G = G1(1,:) + 1i*G1(2,:);      % "

% Sample points
Z = equisphere(0,0,1,M); % M equispaced points in complex unit circle
% Z = randdisk(0,0,1,M); % M random points in complex unit disk
% Z1 = disk_grid(9,1,[0,0],1); % equispaced points in complex unit disk
% Z = Z1(1,:) + 1i*Z1(2,:);    % "
% Z = Z';                      % "


LambdaAll = [];
LambdaCAll = [];
LambdaCAllVF = [];
for m=1:20 % type(m-1,m-1), paper: AAA m <= 16, ratdisk m <= 20

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

rreal = @(x) real(r(x)); % for AAA plot

rfit = @(x) sum(res./(x-pol)); % Vector fitting approximant
rfitreal = @(x) real(rfit(x)); % for vec fitting plot
% optional for plot:
d = abs(fun(0)-rfitreal(0)); % d for Vector fitting
rfitreal = @(x) rfitreal(x) + d; % translation for Vector fitting


figure 
plot(LambdaAll,'o-','Color',[1 0 0])
hold on
leg{1} = 'Condition of Problem';
title_string = strcat('Analytic fct with nearby brance points, $f(z)=\log(1.1 - z)$');
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
title_string = strcat('Analytic fct with nearby brance points, $f(z)=\log(1.1 - z)$');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','Latex')
ylabel('Condition number $\kappa_{AAA}$','Interpreter','Latex')
grid on
hold off

figure
plot(LambdaCAllVF,'o-','Color',[0 1 1]) % Vector fitting
hold on
leg{1} = 'Cauchy matrix VF';
title_string = strcat('Analytic fct with nearby brance points, $f(z)=\log(1.1 - z)$');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
set(gca, 'YScale', 'log')
xlabel('$m$','Interpreter','Latex')
ylabel('Condition number $\kappa_{VF}$','Interpreter','Latex')
grid on
hold off