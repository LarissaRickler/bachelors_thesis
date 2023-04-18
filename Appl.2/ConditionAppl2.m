% Condition of problem and basis of application 2
% 2 Analytic fct with nearby brance points
% log(1.1 - x)

clear all

fun = @(x) log(1.1 -x);
M = 256; % # sample points, paper 256
m = 12; % type(m-1,m-1), paper: AAA m <= 16, ratdisk m <= 20
Z = equisphere(0,0,1,M); % M equispaced points in complex unit circle 
% Z = randdisk(0,0,1,M); % M random points in complex unit disk
% Z1 = disk_grid(9,1,[0,0],1); % equispaced points in complex unit dsk
% Z = Z1(1,:) + 1i*Z1(2,:);
% Z = Z';

% Grid for Lambda
G1 = disk_grid(100,1,[0,0],1);
G = G1(1,:) + 1i*G1(2,:);


[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-13,m);

numL = @(x) sum(abs(w)./abs(x-z)); % numerator of Lebesgue fct
denomL = @(x) abs(sum(w./(x-z))); % denominator of Lebesgue fct
Lambdaf = numL(G)./denomL(G); % Lebesgue fct on points G
Lambda = max(Lambdaf); % "Lebesgue constant" 

C = [];
Zm = setdiff(Z,z);
for k = 1:length(z)
C = [C 1./(Zm-z(k))];
end
LambdaC = cond(C, 2);



% 'Vector fitting' (not proper implemented)
CVF = [];
ZmVF = setdiff(Z,pol);
for k = 1:length(pol)
CVF = [CVF 1./(ZmVF-pol(k))];
end
LambdaCVF = cond(CVF, 2);



rreal = @(x) real(r(x)); % for AAA plot


rfit = @(x) sum(res./(x-pol)); % Vector fitting approximant
rfitreal = @(x) real(rfit(x)); % for vec fitting plot
% optional for plot:
d = abs(fun(0)-rfitreal(0)); % d for Vector fitting
rfitreal = @(x) rfitreal(x) + d; % translation for Vector fitting


%Ratdisk
[rrat,arat,brat,murat,nurat,polesrat,residuesrat] = ratdisk(fun,m-1,m-1,M);
rratreal = @(x) real(rrat(x)); % for plot


figure 
fplot(fun,[-1,1],'-','Color',[1 0 0],'LineWidth',5) % fun
hold on 
fplot(rreal,[-1,1],'-','Color',[0 0 1],'LineWidth',3) % AAA-app
fplot(rratreal,[-1,1],'--','Color',[0 1 0],'LineWidth',3) % ratdisk
fplot(rfitreal,[-1,1],'--','Color',[0 1 1],'LineWidth',3) % VF


leg{1} = 'Function log(1.1 - x)';
leg{2} = 'AAA Approximant';
leg{3} = 'Raddisk Approximant';
leg{4} = 'Vector Fitting Approximant, (not proper implemented)';
title_string = strcat('Analytic fct with nearby brance points');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off

