% Condition of problem and basis of application 3
% 3 Meromorphic fct in the unit disk from boundary values
% tan(beta * x) with beta = 256 in the unit disk

clear all


beta = 256 ; % beta = 4,16,64,256
fun = @(x) tan(beta*x);
M = 1000; % # sample points, paper 1000
m = 200; % type(m-1,m-1), paper: beta/rage(m): 4/(0,14), 16/(0,28), 64/(0,49), 256/(0,62)


% Grid for Lambda
% G1 = disk_grid(100,1,[0,0],1);
% G = G1(1,:) + 1i*G1(2,:);
G = equisphere(0,0,1,3000)';

Z = equisphere(0,0,1,M); % M sample points on unit sphere

[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-13,m);

numL = @(x) sum(abs(w)./abs(x-z)); % numerator of Lebesgue fct
denomL = @(x) abs(sum(w./(x-z))); % denominator of Lebesgue fct
Lambdaf = numL(G)./denomL(G); % Lebesgue fct on points G
Lambda = max(Lambdaf); % "Lebesgue constant" 

C = [];
Zm = setdiff(Z,z);
for k = 1:length(z) % Cauchy matrix
C = [C 1./(Zm-z(k))];
end
LambdaC = cond(C, 2); % condition of Cauchy matrix


% 'Vector fitting'
CVF = [];
ZmVF = setdiff(Z,pol);
for k = 1:length(pol)
CVF = [CVF 1./(ZmVF-pol(k))];
end
LambdaCVF = cond(CVF, 2);


rreal = @(x) real(r(x)); % for AAA plot


rfit = @(x) sum(res./(x-pol)); % Vector fitting approximant (no proper implementation)
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


leg{1} = ['Function tan(' num2str(beta) ' * x)'];
leg{2} = 'AAA Approximant';
leg{3} = 'Raddisk Approximant';
leg{4} = 'Vector Fitting Approximant (no proper implementation)';
title_string = strcat('Meromorphic fct in the unit disk from boundary values');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off



