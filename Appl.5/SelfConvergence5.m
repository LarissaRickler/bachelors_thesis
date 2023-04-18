% Convergence of application 5
% 5 Appr. in connected domains
% Digamma

clear all

M = 2000; % # sample points, paper 2000
fun = @(x) psin(0,x); % Digamma
tol = 1e-14; % 1e-14;  % 1e-13/1e-14

% Sample points
% M random points in complex plane (-5,5)+i(-1,1)
Z = (-5+10*rand(M,1)) + 1i*(-1+2*rand(M,1));                    

% Grid 200x1000 points in complex plane (-5,5)+i(-1,1)
[GX,GY] = meshgrid(linspace(-5,5,1000),linspace(-1,1,200));
G1 = GX + 1i*GY;
G = G1(:)';




AllmaxErrorAAA = [];


for m = 1:20 %max type (m-1,m-1)
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,tol,m); 

errfunAAA = @(x) abs(fun(x)-r(x));

AllmaxErrorAAA = [AllmaxErrorAAA, max(errfunAAA(G))];
end


figure 
plot(AllmaxErrorAAA,'o-','Color',[0 0 1],'LineWidth',2)
hold on
leg{1} = 'AAA Approximant of type(m-1,m-1)';
title_string = strcat('Appr. in connected domains, $f(z) = \psi_0(z)$ with  Digamma func. $\psi_0$');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','LaTex')
ylabel('max error','Interpreter','LaTex')
set(gca, 'YScale', 'log')
grid on
hold off

% levels of error
ReAx = linspace(-10,10,2000);
ImAx = linspace(-3,3,2000);
[RE,IM] = meshgrid(ReAx,ImAx);
errorCont = errfunAAA(RE + 1i*IM);



figure
contour(RE,IM,errorCont,logspace(-12,-1,12)','LineWidth',2);
set(gca,'ColorScale','log')
colormap(parula)
hcb = colorbar;
set(gca, 'clim', [10e-12 10e-1])
hold on 
title_string = strcat('Appr. in connected domains, $f(z) = \psi_0(z)$ with  Digamma func. $\psi_0$');
title(title_string,'Interpreter','LaTex','FontSize',20);
xlabel('Re','Interpreter','LaTex')
ylabel('Im','Interpreter','LaTex')
set(gca,'FontSize',20);
plot(Z,'.','LineWidth',2,'Color','k')
plot(pol,'r*','LineWidth',2)
plot(zer,'g*','LineWidth',2)
leg{1} = 'Max Error';
leg{2} = 'Sample points';
leg{3} = 'Poles of AAA Approximant';
leg{4} = 'Zeros of AAA Approximant';
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
hold off

