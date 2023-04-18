% Convergence of application 6
% 6 Appr. in disconnected domains
% sign(Re(x))

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


AllmaxErrorAAA = [];


for m = 1:100 %max type (m-1,m-1)
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,tol,m); 

errfunAAA = @(x) abs(fun(x)-r(x));

AllmaxErrorAAA = [AllmaxErrorAAA, max(errfunAAA(G))];
end


figure 
plot(AllmaxErrorAAA,'o-','Color',[0 0 1],'LineWidth',2)
hold on
leg{1} = 'AAA Approximant of type(m-1,m-1)';
title_string = strcat(['Appr. in disconnected domains, $f(z) = $sign$($Re$(z))$']);
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
title_string = strcat(['Appr. in disconnected domains, $f(z) = $sign$($Re$(z))$']);
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
