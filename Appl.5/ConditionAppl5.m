% Condition of problem and basis of application 4
% 5 Appr. in connected domains
% 1/bessel

clear all 

M = 2000; % # sample points, paper 2000
m = 13; %max type (m-1,m-1), in paper m = 13
fun = @(x) 1./besselj(0,x); %1/J_0 with J_0 bessel fct

% Grid for Lambda
[GX,GY] = meshgrid(linspace(0,10,1000),linspace(-1,1,200));
G1 = GX + 1i*GY;
G = G1(:)';


LambdaAll = [];
LambdaCAll = [];
LambdaCAllVF = [];
for i=1:100

Z = 10*rand(M,1) + 1i*(-1+2*rand(M,1)); % M random points in complex plane (0,10)+i(-1,1)

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
rfitreal = @(x) real(rfit(x)); % for vec fitting plot


% Ratdisk
[rrat,arat,brat,murat,nurat,polesrat,residuesrat] = ratdisk(fun,m-1,m-1,M);
rratreal = @(x) real(rrat(x)); % for ratdisk plot


figure 
fplot(fun,[0,10],'-','Color',[1 0 0],'LineWidth',5)
hold on 
fplot(rreal,[0,10],'-','Color',[0 0 1],'LineWidth',3) % AAA-App
fplot(rratreal,[0,10],'--','Color',[0 1 0],'LineWidth',3) % ratdisk
fplot(rfitreal,[0,10],'--','Color',[0 1 1],'LineWidth',3) % VF


leg{1} = 'Function 1/bessel';
leg{2} = 'AAA Approximant';
leg{3} = 'Raddisk Approximant';
leg{4} = 'Vector Fitting Approximant';
title_string = strcat('Appr. in connected domains');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off



