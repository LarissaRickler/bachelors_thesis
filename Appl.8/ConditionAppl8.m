% Condition of problem and basis of application 8
% 8 Appr. of exp(x) in (-inf,0] 

clear all;


M = 4000; % # samplepoints
m = 7; % max type (m-1,m-1), paper between 0 and 14
fun = @(x) exp(x);
Z = -logspace(-4,3,M); % M points logarithmically spaced between 10^-4 and 10^3
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-12,m);
% Grid for Lambda
G = -linspace(10^4,-10^(-3),1000000);

numL = @(x) sum(abs(w)./abs(x-z)); % numerator of Lebesgue fct
denomL = @(x) abs(sum(w./(x-z))); % denominator of Lebesgue fct
Lambdaf = numL(G)./denomL(G); % Lebesgue fct on points G
Lambda = max(Lambdaf); % "Lebesgue constant" 

C = [];
Zm = setdiff(Z,z);
for k = 1:length(z)
C = [C; 1./(Zm-z(k))]; % C is different defined
end
C = C';
LambdaC = cond(C, 2);


% 'Vector fitting'
CVF = [];
ZmVF = setdiff(Z,pol);
for k = 1:length(pol)
CVF = [CVF; 1./(ZmVF-pol(k))];
end
CVF = CVF';
LambdaCVF = cond(CVF, 2);


rreal = @(x) real(r(x)); % for AAA plot

rfit = @(x) sum(res./(x-pol)); % Vector fitting approximant (no proper implementation)
rfitreal = @(x) real(rfit(x)); % for vec fitting plot
% optional for plot:
d = abs(fun(0)-rfitreal(0)); % d for Vector fitting
rfitreal = @(x) rfitreal(x) + d; % translation for Vector fitting


% Ratdisk
[rrat,arat,brat,murat,nurat,polesrat,residuesrat] = ratdisk(fun,m-1,m-1,M);
rratreal = @(x) real(rrat(x)); % for ratdisk plot


figure 
fplot(fun,[-100,0],'-','Color',[1 0 0],'LineWidth',5)
hold on 
fplot(rreal,[-100,0],'-','Color',[0 0 1],'LineWidth',3) % AAA-App
fplot(rratreal,[-100,0],'--','Color',[0 1 0],'LineWidth',3) % ratdisk
fplot(rfitreal,[-100,0],'--','Color',[0 1 1],'LineWidth',3) % VF


leg{1} = 'Function';
leg{2} = 'AAA Approximant';
leg{3} = 'Raddisk Approximant';
leg{4} = 'Vector Fitting Approximant (no proper implementation)';
title_string = strcat('Appr. of exp(x) in (-inf,0]');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off





