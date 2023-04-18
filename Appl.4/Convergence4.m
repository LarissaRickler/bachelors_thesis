% Convergence of application 4
% 4 Meromorphic fct in the unit disk from boundary and interior values
% tan(beta * x) with beta = 256 in the unit disk

clear all

M1 = 1000; % # sample points on unit sphere, paper 1000
M2 = 3000; % # sample points in unit disk, paper 3000
M = M1 + M2;
beta = 16; % beta = 4,16,64,256
fun = @(x) tan(beta*x);
tol = 1e-13; % 1e-13/1e-14


Z1 = equisphere(0,0,1,M1); % M1 points in unit sphere
Z2 = randdisk(0,0,1,M2); % M2 points in unit disk
Z = [Z1; Z2]; % sample points 

% Grid 
% G1 = disk_grid(100,1,[0,0],1);
% G = G1(1,:) + 1i*G1(2,:);
G = equisphere(0,0,1,3000);


AllmaxErrorAAA = [];
% AllmaxErrorRat = [];
% AllmaxErrorFit = [];

intpol = [];

for m = 20 % type(m-1,m-1), paper: beta/rage(m): circa 4/(0,14), 16/(0,28), 64/(0,65), 256/(0,200)
% AAA
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,tol,m); 

errfunAAA = @(x) abs(fun(x)-r(x));

AllmaxErrorAAA = [AllmaxErrorAAA, max(errfunAAA(G))];

intpol = [intpol; sum(abs(pol)<1)];
end


figure 
plot(AllmaxErrorAAA,'o-','Color',[0 0 1],'LineWidth',2)
hold on
% plot(AllmaxErrorRat,'o-','Color',[0 1 0])
leg{1} = 'AAA Approximant of type(m-1,m-1)';
% leg{2} = 'Ratdisk Approximant of type(m-1,m-1)';
title_string = strcat(['Meromorphic fct in unit disk from bound and int values, $f(z) = \tan($' num2str(beta) '$z)$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','LaTex')
ylabel('max error','Interpreter','LaTex')
set(gca, 'YScale', 'log')
grid on
hold off



% levels of error
ReAx = linspace(-7.5,7.5,2000);
ImAx = ReAx;
[RE,IM] = meshgrid(ReAx,ImAx);
errorCont = errfunAAA(RE + 1i*IM);

% Unit circle
range = 0:pi/50:2*pi;
xunit = cos(range);
yunit = sin(range);


figure
contour(RE,IM,errorCont,logspace(-12,-1,12)','LineWidth',2);
set(gca,'ColorScale','log')
colormap(parula)
hcb = colorbar;
set(gca, 'clim', [10e-12 10e-1])
hold on 
plot(xunit, yunit,'--','Color', [0 0 0],'LineWidth',2)
title_string = strcat(['Meromorphic fct in unit disk from bound and int values, $f(z) = \tan($' num2str(beta) '$z)$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
xlabel('Re','Interpreter','LaTex')
ylabel('Im','Interpreter','LaTex')
set(gca,'FontSize',20);
plot(pol(2:end),'r*','LineWidth',2)
leg{1} = 'Max Error';
leg{2} = 'Unit Sphere';
leg{3} = 'Poles of AAA Approximant';
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
hold off


