% Condition of problem and basis of application 6
% 6 Appr. in disconnected domains
% sign(Re(x))

clear all

m = 45; % type (m-1,m-1), paper 45-> actually type(43,43), convergence at 51

M = 1000; % # points in each section, paper 1000
Z1 = equisphere(1.5,0,1,M); % M equispaced on unit circle
Z21 = linspace(-2.5,-0.5,M/4 + 1);
Z21 = Z21(1:end-1);
Z22 = linspace(-0.5,-2.5,M/4 + 1);
Z22 = Z22(1:end-1);
Z23 = linspace(1,-1,251);
Z23 = Z23(1:end-1);
Z24 = linspace(-1,1,251);
Z24 = Z24(1:end-1);
Z2 = [ Z21 + 1i, Z22 - 1i, -0.5 + 1i*Z23, -2.5 + 1i*Z24]; % n equispaced around square
Z = [Z1;Z2'];

% % n1 sample points on the unit sphere (equispaced)
% n1 = 1000; % # sample points in circle section
% Z2 = equisphere(1.5,0,1,n1); % n1 equispaced points in unit circle
% % n2^2 sample points in rectangle (equally spaced) 
% n2 = 32;
% [X,Y] = meshgrid(linspace(-2.5,-0.5,n2),linspace(-1,1,n2));
% Z1 = X + 1i*Y;
% Z1 = Z1(:);
% Z = [Z1;Z2];

fun = @(x) sign(real(x));
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-13,m);

% Grid for Lambda
[GX,GY] = meshgrid(linspace(-2.5,-0.5,173),linspace(-1,1,173));
G1 = GX + 1i*GY;
G2 = G1(:)';
G3 = disk_grid(100,1,[1.5,0],10);
G4 = (G3(1,:) + 1i*G3(2,:));
G = [G2, G4];

% Alternative grid for Lambda
[HX,HY] = meshgrid(linspace(-2.5,2.5,1346),linspace(-1,1,173));
H1 = HX + 1i*HY;
H = H1(:)';



numL = @(x) sum(abs(w)./abs(x-z)); % numerator of Lebesgue fct
denomL = @(x) abs(sum(w./(x-z))); % denominator of Lebesgue fct
Lambdaf = numL(G)./denomL(G); % Lebesgue fct on points G
Lambda = max(Lambdaf); % "Lebesgue constant" 


HLambdaf = numL(H)./denomL(H); % Lebesgue fct on points H
HLambda = max(Lambdaf); % "Lebesgue constant" 

C = [];
Zm = setdiff(Z,z);
for k = 1:length(z)
C = [C 1./(Zm-z(k))];
end


% 'Vector fitting'
CVF = [];
ZmVF = setdiff(Z,pol);
for k = 1:length(pol)
CVF = [CVF 1./(ZmVF-pol(k))];
end
LambdaCVF = cond(CVF, 2);
LambdaC = cond(C, 2);


rreal = @(x) real(r(x)); % for AAA plot

rfit = @(x) sum(res./(x-pol)); % Vector fitting approximant
rfitreal1 = @(x) real(rfit(x)); % for vec fitting plot
% optional for plot:
d = abs(rfitreal1(1)-fun(1)); % d for Vector fitting
rfitreal = @(x) rfitreal1(x) - d; % translation for Vector fitting


% % Ratdisk, not in the unit disk
% [rrat,arat,brat,murat,nurat,polesrat,residuesrat] = ratdisk(fun,m-1,m-1,2*M+1);
% rratreal = @(x) real(rrat(x)); % for ratdisk plot


figure 
fplot(fun,[-2.5,2.5],'-','Color',[1 0 0],'LineWidth',5)
hold on 
fplot(rreal,[-2.5,2.5],'-','Color',[0 0 1],'LineWidth',3) % AAA-App
% fplot(rratreal,[-2.5,2.5],'--','Color',[0 1 0],'LineWidth',3) % ratdisk
fplot(rfitreal,[-2.5,2.5],'--','Color',[0 1 1],'LineWidth',3) % VF


leg{1} = 'Function sign(Re(x))';
leg{2} = 'AAA Approximant';
% leg{3} = 'Raddisk Approximant';
leg{4} = 'Vector Fitting Approximant';
title_string = strcat('Appr. in disconnected domains');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off





