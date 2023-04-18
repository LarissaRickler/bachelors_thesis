% Condition of problem and basis of application 7
% 7 Approximation of |x| on [-1,1]
% Way two

clear all

M = 100000; % # samplepoints in [0,1] (half as much as in 'way one')
m = 3; % max type (m-1,m-1), (circa half of 'way one' for equal good approx.) 
fun1 = @(x) sqrt(x);
Z = linspace(0,1,M); % M points equispaced in [0,1]   
[r1,pol,res,zer,z,f,w,errvec] = aaa(fun1(Z),Z,1e-13,m);
%convert
fun = @(x) fun1(x.^2);
r = @(x) r1(x.^2);

% Grid for Lambda 
% G = [-Z,Z];
G = [Z];


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

rfit = @(x) sum(res./(x.^2 - pol)); % Vector fitting approximant (no proper implementation)
rfitreal = @(x) real(rfit(x)); % for vec fitting plot
% optional for plot:
d = abs(fun(0)-rfitreal(0)); % d for Vector fitting
rfitreal = @(x) rfitreal(x) + d; % translation for Vector fitting


% Ratdisk
[rrat,arat,brat,murat,nurat,polesrat,residuesrat] = ratdisk(fun1,m-1,m-1,M);
rratreal = @(x) real(rrat(x.^2)); % for ratdisk plot


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





