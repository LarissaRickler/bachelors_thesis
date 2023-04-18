% Condition of problem and basis of application 4
% 4 Meromorphic fct in the unit disk from bound and int values
% tan(beta * x) with beta = 256 in the unit disk from boundary and
%  interior values
% 100 trys to compare the Conditions, since Z is random

clear all

M1 = 1000; % # sample points on unit sphere, paper 1000
M2 = 3000; % # sample points in unit disk, paper 3000
M = M1 + M2;
m = 200; % type(m-1,m-1), paper: beta/rage(m): circa 4/(0,14), 16/(0,28), 64/(0,65), 256/(0,200)
Z1 = equisphere(0,0,1,M1); % M1 points in unit sphere
beta = 4; % beta in paper 4,16,64,256 
fun = @(x) tan(beta*x);

% Grid for Lambda
G1 = disk_grid(100,1,[0,0],1);
G = G1(1,:) + 1i*G1(2,:);

LambdaAll = [];
LambdaCAll = [];
LambdaCAllVF = [];
for i=1:100

Z2 = randdisk(0,0,1,M2); % M2 points in unit disk
Z = [Z1; Z2]; % sample points 

[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-13,m);

numL = @(x) sum(abs(w)./abs(x-z)); % numerator of Lebesgue fct
denomL = @(x) abs(sum(w./(x-z))); % denominator of Lebesgue fct
Lambdaf = numL(G)./denomL(G); % Lebesgue fct on points G
Lambda = max(Lambdaf); % "Lebesgue constant" 
LambdaAll = [LambdaAll, Lambda];
C = [];
Zm = setdiff(Z,z);
for k = 1:length(z)
C = [C 1./(Zm-z(k))];
end
LambdaC = cond(C, 2);
LambdaCAll = [LambdaCAll, LambdaC];


% 'Vector fitting'
CVF = [];
ZmVF = setdiff(Z,pol);
for k = 1:length(pol)
CVF = [CVF 1./(ZmVF-pol(k))];
end
LambdaCVF = cond(CVF, 2);
LambdaCAllVF = [LambdaCAllVF, LambdaCVF];
end
disp('Done')


rreal = @(x) real(r(x)); % for AAA plot


rfit = @(x) sum(res./(x-pol)); % Vector fitting approximant
rfitreal1 = @(x) real(rfit(x)); % for vec fitting plot
% optional for plot:
d = abs(fun(0)-rfitreal1(0)); % d for Vector fitting
rfitreal = @(x) rfitreal1(x) + d; % translation for Vector fitting


%Ratdisk
[rrat,arat,brat,murat,nurat,polesrat,residuesrat] = ratdisk(fun,m-1,m-1,M1);
rratreal = @(x) real(rrat(x)); % for plot


figure 
fplot(fun,[-1,1],'-','Color',[1 0 0],'LineWidth',5) % fun
hold on 
fplot(rreal,[-1,1],'-','Color',[0 0 1],'LineWidth',3) % AAA-app
fplot(rratreal,[-1,1],'--','Color',[0 1 0],'LineWidth',3) % ratdisk
fplot(rfitreal,[-1,1],'--','Color',[0 1 1],'LineWidth',3) % VF


leg{1} = 'Function tan(beta * x)';
leg{2} = 'AAA Approximant';
leg{3} = 'Raddisk Approximant';
leg{4} = 'Vector Fitting Approximant';
title_string = strcat('Meromorphic fct in the unit disk from bound and int values');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off


