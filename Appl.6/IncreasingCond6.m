% Increasing depending on type
% Condition of problem and basis of application 6
% 6 Appr. in disconnected domains
% sign(real(x))

clear all


fun = @(x) sign(real(x));
tol = 1e-14; % tolerance 1e-13/1e-14


% Sample points
M = 1000; % # points in each section, paper 1000
Z1 = equisphere(1.5,0,1,M); % M equispaced on unit circle
Z21 = linspace(-2.5,-0.5,M/4 + 1);
Z21 = Z21(1:end-1);
Z22 = linspace(-0.5,-2.5,M/4 + 1);
Z22 = Z22(1:end-1);
Z23 = linspace(1,-1,M/4 + 1);
Z23 = Z23(1:end-1);
Z24 = linspace(-1,1,M/4 + 1);
Z24 = Z24(1:end-1);
Z2 = [ Z21 + 1i, Z22 - 1i, -0.5 + 1i*Z23, -2.5 + 1i*Z24]; % n equispaced around square
Z = [Z1;Z2'];
                

% Grid 173x173 points in complex plane (-2.5,-0.5)+i(-1,1)
% and 31757 equispaced points in disk with radius 1 and center [1.5,0]
[GX,GY] = meshgrid(linspace(-2.5,-0.5,173),linspace(-1,1,173));
G1 = GX + 1i*GY;
G2 = G1(:)';
G3 = disk_grid(100,1,[1.5,0],10);
G4 = (G3(1,:) + 1i*G3(2,:));
G = [G2, G4];



LambdaAll = [];
LambdaCAll = [];
LambdaCAllVF = [];
for m=1:50 %max type (m-1,m-1), in paper m = 15

[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,tol,m);

numL = @(x) sum(abs(w)./abs(x-z)); % numerator of Lebesgue fct
denomL = @(x) abs(sum(w./(x-z))); % denominator of Lebesgue fct
Lambdaf = numL(G)./denomL(G); % Lebesgue fct on points G
Lambda = max(Lambdaf); % "Lebesgue constant" 
LambdaAll = [LambdaAll, Lambda];

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
end


figure 
plot(LambdaAll,'o-','Color',[1 0 0]) % Lambda
hold on
leg{1} = 'Condition of Problem';
title_string = strcat(['Appr. in disconnected domains, $f(z) = $sign$($Re$(z))$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
set(gca, 'YScale', 'log')
xlabel('$m$','Interpreter','Latex')
ylabel('Condition number $\Lambda^{(m)}$','Interpreter','Latex')
grid on
hold off

figure 
plot(LambdaCAll,'o-','Color',[0 0 1]) % AAA-App
hold on
leg{1} = 'Cauchy matrix AAA';
title_string = strcat(['Appr. in disconnected domains, $f(z) = $sign$($Re$(z))$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','Latex')
ylabel('Condition number $\kappa_{AAA}$','Interpreter','Latex')
grid on
hold off

figure
plot(LambdaCAllVF,'o-','Color',[0 1 1]) % VF
hold on
leg{1} = 'Cauchy matrix VF';
title_string = strcat(['Appr. in disconnected domains, $f(z) = $sign$($Re$(z))$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
set(gca, 'YScale', 'log')
xlabel('$m$','Interpreter','Latex')
ylabel('Condition number $\kappa_{VF}$','Interpreter','Latex')
grid on
hold off