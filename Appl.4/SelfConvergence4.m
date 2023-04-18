% Convergence of application 4
% 4 Meromorphic fct in the unit disk from boundary values
% f(z) = 100.*sec(beta.*x).*csc(beta.*x), beta = 4,16 

clear all

beta1 = 4; % beta = 4,16 -> 5 and 21 poles
fun1 = @(x) (100).*sec(beta1.*x).*csc(beta1.*x); 

beta2 = 16; % beta = 4,16 -> 5 and 21 poles
fun2 = @(x) (100).*sec(beta2.*x).*csc(beta2.*x);

random = 0; % set true for random sample set
NoSpuriousPoles = 0; % set true for working space with no spurious poles
SpuriousPoles = 1; % set true for working space with spurious poles

% Grid 
G1 = disk_grid(100,1,[0,0],1); % equispaced points in complex unit disk
G = G1(1,:) + 1i*G1(2,:);      % "
GB = equisphere(0,0,1,3100)';% equispaced points in complex unit sphere

tol = 1e-14; % 1e-13/1e-14 tolerance


if(random)
    M1 = 501; % # sample points on unit sphere
    M21 = 3000; % # sample points in unit disk for beta = 4
    M22 = 3000; % # sample points in unit disk for beta = 16
    Mb4 = M1 + M21; % # sample points for beta = 4
    Mb16 = M1 + M22; % # sample points for beta = 16

    % Sample points
    Z1 = equisphere(0,0,1,M1); % M1 points in unit sphere
    Z21 = randdisk(0,0,1,M21); % M2 points in unit disk
    Z22 = randdisk(0,0,1,M22); % M2 points in unit disk
    Zb4 = [Z1; Z21]; % sample points for beta = 4
    Zb16 = [Z1; Z22]; % sample points for beta = 16
end 

if(NoSpuriousPoles)
    load('WorkspaceWithoutSpuriousPoles')
end

if(SpuriousPoles)
    load('WorkspaceWithSpuriousPoles')
end




AllmaxErrorAAA1 = [];
AllmaxErrorAAA2 = [];

intpol1 = [];
intpol2 = [];

for m = 1:50 % type(m-1,m-1) at m = 10 winding number -5 for beta = 4
% % AAA
[r1,pol1,res1,zer1,z1,f1,w1,errvec1] = aaa(fun1(Zb4),Zb4,tol,m); 
[r2,pol2,res2,zer2,z2,f2,w2,errvec2] = aaa(fun2(Zb16),Zb16,tol,m);



errfunAAA1 = @(x) abs(fun1(x)-r1(x));
errfunAAA2 = @(x) abs(fun2(x)-r2(x));

AllmaxErrorAAA1 = [AllmaxErrorAAA1, max(errfunAAA1(GB))];
AllmaxErrorAAA2 = [AllmaxErrorAAA2, max(errfunAAA2(GB))];

intpol1 = [intpol1; sum(abs(pol1)<1)];
intpol2 = [intpol2; sum(abs(pol2)<1)];
end





figure 
subplot(1,2,1)
plot(AllmaxErrorAAA1,'o-','Color',[0 0 1],'LineWidth',2)
hold on
leg{1} = 'AAA Approximant of type(m-1,m-1)';
title_string = strcat(['Meromorphic fct in the unit disk from boundary and int values, $f(z) = 100\sec($' num2str(beta1) '$z)\csc($' num2str(beta1) '$z)$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','LaTex')
ylabel('max error','Interpreter','LaTex')
set(gca, 'YScale', 'log')
grid on

subplot(1,2,2)
plot(AllmaxErrorAAA2,'o-','Color',[0 0 1],'LineWidth',2)
hold on
leg{1} = 'AAA Approximant of type(m-1,m-1)';
title_string = strcat(['Meromorphic fct in the unit disk from boundary and int values, $f(z) = 100\sec($' num2str(beta2) '$z)\csc($' num2str(beta2) '$z)$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','LaTex')
ylabel('max error','Interpreter','LaTex')
set(gca, 'YScale', 'log')
grid on
hold off



mplot1 = 19; % type (mplot1-1,mplot1-1) for cotour plot beta = 4
mplot2 = 38; % type (mplot2-1,mplot2-1) for cotour plot beta = 16

% AAA
[r1,pol1,res1,zer1,z1,f1,w1,errvec1] = aaa(fun1(Zb4),Zb4,tol,mplot1); 
[r2,pol2,res2,zer2,z2,f2,w2,errvec2] = aaa(fun2(Zb16),Zb16,tol,mplot2);


errfunAAA1 = @(x) abs(fun1(x)-r1(x));

errfunAAA2 = @(x) abs(fun2(x)-r2(x));



% levels of error
ReAx = linspace(-3,3,2000);
ImAx = ReAx;
[RE,IM] = meshgrid(ReAx,ImAx);
errorCont1 = errfunAAA1(RE + 1i*IM);
errorCont2 = errfunAAA2(RE + 1i*IM);



% Unit circle
range = 0:pi/50:2*pi;
xunit = cos(range);
yunit = sin(range);


figure
subplot(1,2,1)
contour(RE,IM,errorCont1,logspace(-12,-1,12)','LineWidth',2);
set(gca,'ColorScale','log')
colormap(parula)
hcb = colorbar;
set(gca, 'clim', [10e-12 10e-1])
hold on 
plot(xunit, yunit,'--','Color', [0 0 0],'LineWidth',2)
title_string = strcat(['Meromorphic fct in the unit disk from boundary and int values, $f(z) = 100\sec($' num2str(beta1) '$z)\csc($' num2str(beta1) '$z)$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
xlabel('Re','Interpreter','LaTex')
ylabel('Im','Interpreter','LaTex')
set(gca,'FontSize',20);
plot(pol1,'r*','LineWidth',2)
leg{1} = ['Max Error of AAA Approximant of type $($' num2str(mplot1-1) '$,$' num2str(mplot1-1) '$)$'];
leg{2} = 'Unit Sphere';
leg{3} = 'Poles of AAA Approximant';
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
hold off

subplot(1,2,2)
contour(RE,IM,errorCont2,logspace(-12,-1,12)','LineWidth',2);
set(gca,'ColorScale','log')
colormap(parula)
hcb = colorbar;
set(gca, 'clim', [10e-12 10e-1])
hold on 
plot(xunit, yunit,'--','Color', [0 0 0],'LineWidth',2)
title_string = strcat(['Meromorphic fct in the unit disk from boundary and int values, $f(z) = 100\sec($' num2str(beta2) '$z)\csc($' num2str(beta2) '$z)$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
xlabel('Re','Interpreter','LaTex')
ylabel('Im','Interpreter','LaTex')
set(gca,'FontSize',20);
plot(pol2,'r*','LineWidth',2)
leg{1} = ['Max Error of AAA Approximant of type $($' num2str(mplot2-1) '$,$' num2str(mplot2-1) '$)$']';
leg{2} = 'Unit Sphere';
leg{3} = 'Poles of AAA Approximant';
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
hold off

