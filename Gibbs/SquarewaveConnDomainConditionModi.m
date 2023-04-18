% Increasing depending on type
% Condition of problem and basis of approximation discontinuous functions
% Modified
% Square wave

clear all

fun = @(x) sign(sin((2.*pi.*real(x))./2));
tol = 1e-13; % tolerance 1e-13/1e-14


S = @(x) sGibbs([1 2 3 4],[2 2 2 2],10,x);


M = 5000; % Sample points per domain
Z = linspace(0,5,M);
Z = Z(2:end);
SZ = S(Z);

G = linspace(0,5,5000);
G = setdiff(G,[0,1,2,3,4])';



LambdaAll = [];
LambdaCAll = [];
LambdaCAllVF = [];
for m= 1:100%max type (m-1,m-1), in paper m = 15

[r1,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),SZ,tol,m);
r = @(x) r1(S(x));

numL = @(x) sum(abs(w)./abs(x-z)); % numerator of Lebesgue fct
denomL = @(x) abs(sum(w./(x-z))); % denominator of Lebesgue fct
Lambdaf = numL(S(G)')./denomL(S(G)'); % Lebesgue fct on points G
Lambda = max(Lambdaf); % "Lebesgue constant" 
LambdaAll = [LambdaAll, Lambda];

C = [];
Zm = setdiff(SZ,z);
for k = 1:length(z) % Cauchy matrix AAA
C = [C 1./(Zm-z(k))];
end

LambdaC = cond(C, 2);
LambdaCAll = [LambdaCAll, LambdaC];


% 'Vector fitting'
CVF = [];
ZmVF = setdiff(SZ,pol);
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
title_string = strcat(['Appr. discontinuous functions, Modified, $f(z) = $ sign$(\sin(\frac{2\pi}{2}$Re$(z)))$']);
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
title_string = strcat(['Appr. discontinuous functions, Modified, $f(z) = $ sign$(\sin(\frac{2\pi}{2}$Re$(z)))$']);
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
title_string = strcat(['Appr. discontinuous functions, Modified, $f(z) = $ sign$(\sin(\frac{2\pi}{2}$Re$(z)))$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
set(gca, 'YScale', 'log')
xlabel('$m$','Interpreter','Latex')
ylabel('Condition number $\kappa_{VF}$','Interpreter','Latex')
grid on
hold off