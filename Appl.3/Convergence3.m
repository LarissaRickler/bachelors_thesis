% Convergence of application 3
% 3 Meromorphic fct in the unit disk from boundary values
% tan(beta * x) in the unit disk

clear all

beta = 64; % beta = 4,16,64,256
fun = @(x) tan(beta*x);
M = 1000; % # sample points, paper 1000
tol = 1e-13; % 1e-13/1e-14 tolerance

% Sample points
Z = equisphere(0,0,1,M); % M equispaced points in complex unit circle  

% Grid 
% G1 = disk_grid(100,1,[0,0],1); % equispaced points in complex unit disk
% G = G1(1,:) + 1i*G1(2,:);      % "
G = equisphere(0,0,1,3100); % equispaced points in complex unit sphere




AllmaxErrorAAA = [];
AllmaxErrorRat = [];
AllmaxErrorPoly = [];

for m=1:70 % type(m-1,m-1), paper: beta/rage(m): 4/(0,14), 16/(0,28), 64/(0,49), 256/(0,62)
% AAA
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,tol,m); 


% Ratdisk
[rrat,arat,brat,murat,nurat,polesrat,residuesrat] = ratdisk(fun,m-1,m-1,M,tol);


% Polynom interpolation
pc = polyfit(Z,fun(Z),2*(m-1));
p = @(x) polyval(pc,x);


errfunAAA = @(x) abs(fun(x)-r(x));
errfunrat  = @(x) abs(fun(x)-rrat(x));
errfunpoly = @(x) abs(fun(x)-p(x));


AllmaxErrorAAA = [AllmaxErrorAAA, max(errfunAAA(G))];
AllmaxErrorRat = [AllmaxErrorRat, max(errfunrat(G))];
AllmaxErrorPoly = [AllmaxErrorPoly, max(errfunpoly(G))];
end


figure 
plot(AllmaxErrorAAA,'o-','Color',[0 0 1],'LineWidth',2)
hold on
plot(AllmaxErrorRat,'o-','Color',[0 1 0],'LineWidth',2)
plot(AllmaxErrorPoly,'o-','Color',[1 0 1],'LineWidth',2)
leg{1} = 'AAA Approximant of type(m-1,m-1)';
leg{2} = 'Ratdisk Approximant of type(m-1,m-1)';
leg{3} = 'Polynomial of degree 2(m-1)';
title_string = strcat(['Meromorphic fct in the unit disk from boundary values, $f(z) = \tan($' num2str(beta) '$z)$']);
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
title_string = strcat(['Meromorphic fct in the unit disk from boundary values, $f(z) = \tan($' num2str(beta) '$z)$']);
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


