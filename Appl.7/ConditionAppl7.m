% Condition of problem and basis of application 7
% 7 Approximation of |x| on [-1,1]
% Way one

clear all

M = 200000; % # samplepoints in [-1,1], paper 200000
m = 5; % max type (m-1,m-1), paper between 0 and 26
fun = @(x) abs(x);
Z = linspace(-1,1,M); % M equi spaced points in [-1,1]
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-13,m);
% Grid for Lambda
G = Z;


numL = @(x) sum(abs(w)./abs(x-z)); % numerator of Lebesgue fct
denomL = @(x) abs(sum(w./(x-z))); % denominator of Lebesgue fct
Lambdaf = numL(G)./denomL(G); % Lebesgue fct on points G
Lambda = max(Lambdaf); % "Lebesgue constant" 

C = [];
Zm = setdiff(Z,z);
for k = 1:length(z)
C = [C; 1./(Zm-z(k))]; % C is different defined
end
LambdaC = cond(C, 2);


% 'Vector fitting'
CVF = [];
ZmVF = setdiff(Z,pol);
for k = 1:length(pol)
CVF = [CVF; 1./(ZmVF-pol(k))];
end
LambdaCVF = cond(CVF, 2);


rreal = @(x) real(r(x)); % for AAA plot

rfit = @(x) sum(res./(x-pol)); % Vector fitting approximant (no proper implementation)
rfitreal = @(x) real(rfit(x)); % for vec fitting plot
% optional for plot:
d = abs(fun(0)-rfitreal(0)); % d for Vector fitting
rfitreal = @(x) rfitreal(x) + d; % translation for Vector fitting


% Ratdisk
[rrat,arat,brat,murat,nurat,polesrat,residuesrat] = ratdisk(fun,m-1,m-1,M);
rratreal = @(x) real(rrat(2.*x)); % for ratdisk plot


figure 
fplot(fun,[-1,1],'-','Color',[1 0 0],'LineWidth',5)
hold on 
fplot(rreal,[-1,1],'-','Color',[0 0 1],'LineWidth',3) % AAA-App
fplot(rratreal,[-1,1],'--','Color',[0 1 0],'LineWidth',3) % ratdisk
fplot(rfitreal,[-1,1],'--','Color',[0 1 1],'LineWidth',3) % VF


leg{1} = 'Function';
leg{2} = 'AAA Approximant';
leg{3} = 'Raddisk Approximant';
leg{4} = 'Vector Fitting Approximant (no proper implementation)';
title_string = strcat('Approximation of |x| on [-1,1]');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off





