% Condition of problem and basis of application 9
% 9 Clamped beam model form Chahlaoui and Van Dooren
% fun(z) = c'(zI-A)^(-1)b (a rat. fun. of type(348,348) 

clear all

load('ClampedBeam.mat')

M = 500; % % 2*M sample points, paper n = 500
m = 3; % type (m-1,m-1), paper between 0 and 46
Z1 = 1i*logspace(-2,2,M); % M logarithmically spaced points from 10?2 i to 102 i
Z2 = -Z1; % complex conjugates
dim = 348;
Z = [Z1 Z2];
fun = @(x) cAb(C,A,B,x);
funZ = fun(Z);

[r,pol,res,zer,z,f,w,errvec] = aaa(funZ,Z,1e-13,m);

% Grid for Lambda
G = 1i*linspace(-100,100,20000);
funG = fun(G);


numL = @(x) sum(abs(w)./abs(x-z)); % numerator of Lebesgue fct
denomL = @(x) abs(sum(w./(x-z))); % denominator of Lebesgue fct
Lambdaf = numL(G)./denomL(G); % Lebesgue fct on points G
Lambda = max(Lambdaf); % "Lebesgue constant" 

CAAA = [];
Zm = setdiff(Z,z);
for k = 1:length(z)
CAAA = [CAAA; 1./(Zm-z(k))]; % C is different defined
end
CAAA = CAAA';
LambdaC = cond(CAAA, 2);


% 'Vector fitting'
CVF = [];
ZmVF = setdiff(Z,pol);
for k = 1:length(pol)
CVF = [CVF; 1./(ZmVF-pol(k))];
end
CVF = CVF';
LambdaCVF = cond(CVF, 2);


rreal = @(x) real(r(x)); % for AAA plot

rfit = @(x) sum(res./(x-pol)); % Vector fitting approximant
rfitreal = @(x) real(rfit(x)); % for vec fitting plot
% optional for plot:
d = abs(fun(0)-rfitreal(0)); % d for Vector fitting
rfit = @(x) rfit(x) + d;
rfitreal = @(x) rfitreal(x) + d; % translation for Vector fitting


figure 
plot(imag(G),abs(fun(G)),'-','Color',[1 0 0],'LineWidth',5)
hold on 
plot(imag(G),abs(r(G)),'-','Color',[0 0 1],'LineWidth',3) % AAA-App
plot(imag(G),abs(rfit(G)),'--','Color',[0 1 1],'LineWidth',3) % VF
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')


leg{1} = 'Function';
leg{2} = 'AAA Approximant';
leg{3} = 'Raddisk Approximant';
leg{4} = 'Vector Fitting Approximant';
title_string = strcat('Clamped beam model form Chahlaoui and Van Dooren');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off




