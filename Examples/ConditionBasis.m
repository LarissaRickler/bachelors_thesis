% Condition of basis

% % Example 1
% clear all
% 
% M = 2000; % # sample points
% fun = @(x) sin(x).*x.^3; 
% tol = 1e-14; % 1e-13/1e-14
% 
% % Grid for Lambda
% [GX,GY] = meshgrid(linspace(0,10,1000),linspace(-1,1,200));
% G1 = GX + 1i*GY;
% G = G1(:)';
% 
% Z = 10*rand(M,1) + 1i*(-1+2*rand(M,1)); % M random points in complex plane (0,10)+i(-1,1)


% % Example 2
% clear all
% 
% M = 2000; % # sample points, paper 2000
% fun = @(x) psin(0,x); % Digamma
% tol = 1e-14; % 1e-14;  % 1e-13/1e-14
% 
% % Sample points
% M random points in complex plane (-5,5)+i(-1,1)
% Z = (-5+10*rand(M,1)) + 1i*(-1+2*rand(M,1));                    
% 
% % Grid 200x1000 points in complex plane (-5,5)+i(-1,1)
% [GX,GY] = meshgrid(linspace(-5,5,1000),linspace(-1,1,200));
% G1 = GX + 1i*GY;
% G = G1(:)';


LambdaCAll = [];
LambdaCAllVF = [];
LambdaVan = [];
for m=1:16 % type (m-1,m-1) 

[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,tol,m);

C = [];
Zm = setdiff(Z,z);
for k = 1:length(z) % Cauchy matrix AAA
C = [C 1./(Zm-z(k))];
end
LambdaC = cond(C, 2);
LambdaCAll = [LambdaCAll, LambdaC];


% 'Vector fitting'
CVF = [];
ZmVF = setdiff(Z,pol);
for k = 1:length(pol) % Cauchy matrix VF
CVF = [CVF 1./(ZmVF-pol(k))];
end
LambdaCVF = cond(CVF, 2);
LambdaCAllVF = [LambdaCAllVF, LambdaCVF];


Van = fliplr(vander(Z));
VanC = cond(Van(:,1:m),2);
LambdaVan = [LambdaVan, VanC];
end
LambdaCAll = LambdaCAll(1:end-1);
LambdaCAllVF = LambdaCAllVF(2:end);
LambdaVan = LambdaVan(1:end-1);


figure 
plot(LambdaCAll,'o-','Color',[0 0 1],'LineWidth',2) % AAA-App
hold on
plot(LambdaVan,'o-','Color',[0 1 0],'LineWidth',2) % Ratdisk
plot(LambdaCAllVF,'o-','Color',[0 1 1],'LineWidth',2) % VF
leg{1} = 'Barycentric basis';
leg{2} = 'Monomial basis';
leg{3} = 'Partial fraction basis';
legend(leg,'Interpreter','LaTex','FontSize',14,'Location','EastOutside');
set(gca,'FontSize',14);
set(gca, 'YScale', 'log')
set(gca,'xtick',[0:2:16])
xlabel('$m$','Interpreter','Latex')
ylabel('Condition number','Interpreter','Latex')
grid on
hold off
