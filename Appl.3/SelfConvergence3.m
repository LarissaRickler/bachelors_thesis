% Convergence of application 3
% 3 Meromorphic fct in the unit disk from boundary values
% f(z) = 100.*sec(beta.*x).*csc(beta.*x), beta = 4,16 

clear all



M = 501; % # sample points

beta1 = 4; % beta = 4,16 -> 5 and 21 poles
fun1 = @(x) (100).*sec(beta1.*x).*csc(beta1.*x); 

beta2 = 16; % beta = 4,16 -> 5 and 21 poles
fun2 = @(x) (100).*sec(beta2.*x).*csc(beta2.*x);

tol = 1e-14; % 1e-13/1e-14 tolerance

% Sample points
Z = equisphere(0,0,1,M); % M equispaced points in complex unit circle  

% Grid 
% G1 = disk_grid(100,1,[0,0],1); % equispaced points in complex unit disk
% G = G1(1,:) + 1i*G1(2,:);      % "
G = equisphere(0,0,1,3100)';% equispaced points in complex unit sphere





AllmaxErrorAAA1 = [];
AllmaxErrorRat1 = [];
AllmaxErrorAAA2 = [];
AllmaxErrorRat2 = [];

intpol1 = [];
intpol2 = [];
intpolrat1 = [];
intpolrat2 = [];

for m = 1:50 % type(m-1,m-1) at m = 10 winding number -5 for beta = 4
% % AAA
[r1,pol1,res1,zer1,z1,f1,w1,errvec1] = aaa(fun1(Z),Z,tol,m); 
[r2,pol2,res2,zer2,z2,f2,w2,errvec2] = aaa(fun2(Z),Z,tol,m);

% % Ratdisk
[rrat1,arat1,brat1,murat1,nurat1,polesrat1,residuesrat1] = ratdisk(fun1,m-1,m-1,M,tol);
[rrat2,arat2,brat2,murat2,nurat2,polesrat2,residuesrat2] = ratdisk(fun2,m-1,m-1,M,tol);

errfunAAA1 = @(x) abs(fun1(x)-r1(x));
errfunrat1  = @(x) abs(fun1(x)-rrat1(x));
errfunAAA2 = @(x) abs(fun2(x)-r2(x));
errfunrat2  = @(x) abs(fun2(x)-rrat2(x));

AllmaxErrorAAA1 = [AllmaxErrorAAA1, max(errfunAAA1(G))];
AllmaxErrorRat1 = [AllmaxErrorRat1, max(errfunrat1(G))];
AllmaxErrorAAA2 = [AllmaxErrorAAA2, max(errfunAAA2(G))];
AllmaxErrorRat2 = [AllmaxErrorRat2, max(errfunrat2(G))];

intpol1 = [intpol1; sum(abs(pol1)<1)];
intpol2 = [intpol2; sum(abs(pol2)<1)];

intpolrat1 = [intpolrat1, sum(abs(polesrat1)<1)];
intpolrat2 = [intpolrat2, sum(abs(polesrat2)<1)];
end




figure 
subplot(1,2,1)
plot(AllmaxErrorAAA1,'o-','Color',[0 0 1],'LineWidth',2)
hold on
plot(AllmaxErrorRat1,'o-','Color',[0 1 0],'LineWidth',2)
leg{1} = 'AAA Approximant of type(m-1,m-1)';
leg{2} = 'Ratdisk Approximant of type(m-1,m-1)';
title_string = strcat(['Meromorphic fct in the unit disk from boundary values, $f(z) = 100\sec($' num2str(beta1) '$z)\csc($' num2str(beta1) '$z)$']);
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
plot(AllmaxErrorRat2,'o-','Color',[0 1 0],'LineWidth',2)
leg{1} = 'AAA Approximant of type(m-1,m-1)';
leg{2} = 'Ratdisk Approximant of type(m-1,m-1)';
title_string = strcat(['Meromorphic fct in the unit disk from boundary values, $f(z) = 100\sec($' num2str(beta2) '$z)\csc($' num2str(beta2) '$z)$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','LaTex')
ylabel('max error','Interpreter','LaTex')
set(gca, 'YScale', 'log')
grid on
hold off



mplot1 = 20; % type (mplot1-1,mplot1-1) for cotour plot beta = 4
mplot2 = 40; % type (mplot2-1,mplot2-1) for cotour plot beta = 16

% AAA
[r1,pol1,res1,zer1,z1,f1,w1,errvec1] = aaa(fun1(Z),Z,tol,mplot1); 
[r2,pol2,res2,zer2,z2,f2,w2,errvec2] = aaa(fun2(Z),Z,tol,mplot2);

% Ratdisk
[rrat1,arat1,brat1,murat1,nurat1,polesrat1,residuesrat1] = ratdisk(fun1,mplot1-1,mplot1-1,M,tol);
[rrat2,arat2,brat2,murat2,nurat2,polesrat2,residuesrat2] = ratdisk(fun2,mplot2-1,mplot2-1,M,tol);

errfunAAA1 = @(x) abs(fun1(x)-r1(x));
errfunrat1  = @(x) abs(fun1(x)-rrat1(x));
errfunAAA2 = @(x) abs(fun2(x)-r2(x));
errfunrat2  = @(x) abs(fun2(x)-rrat2(x));


% levels of error
ReAx = linspace(-3,3,2000);
ImAx = ReAx;
[RE,IM] = meshgrid(ReAx,ImAx);
errorCont1 = errfunAAA1(RE + 1i*IM);
errorCont2 = errfunAAA2(RE + 1i*IM);
errorContrat1 = errfunrat1(RE + 1i*IM);
errorContrat2 = errfunrat2(RE + 1i*IM);


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
title_string = strcat(['Meromorphic fct in the unit disk from boundary values, $f(z) = 100\sec($' num2str(beta1) '$z)\csc($' num2str(beta1) '$z)$']);
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
title_string = strcat(['Meromorphic fct in the unit disk from boundary values, $f(z) = 100\sec($' num2str(beta2) '$z)\csc($' num2str(beta2) '$z)$']);
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

polesrat1 = complex(polesrat1);

figure
subplot(1,2,1)
contour(RE,IM,errorContrat1,logspace(-12,-1,12)','LineWidth',2);
set(gca,'ColorScale','log')
colormap(parula)
hcb = colorbar;
set(gca, 'clim', [10e-12 10e-1])
hold on 
plot(xunit, yunit,'--','Color', [0 0 0],'LineWidth',2)
title_string = strcat(['Meromorphic fct in the unit disk from boundary values, $f(z) = 100\sec($' num2str(beta1) '$z)\csc($' num2str(beta1) '$z)$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
xlabel('Re','Interpreter','LaTex')
ylabel('Im','Interpreter','LaTex')
set(gca,'FontSize',20);
plot(polesrat1,'r*','LineWidth',2)
leg{1} = ['Max Error of Ratdisk Approximant of type $($' num2str(mplot1-1) '$,$' num2str(mplot1-1) '$)$']';
leg{2} = 'Unit Sphere';
leg{3} = 'Poles of Ratdisk Approximant';
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
hold off

subplot(1,2,2)
contour(RE,IM,errorContrat2,logspace(-12,-1,12)','LineWidth',2);
set(gca,'ColorScale','log')
colormap(parula)
hcb = colorbar;
set(gca, 'clim', [10e-12 10e-1])
hold on 
plot(xunit, yunit,'--','Color', [0 0 0],'LineWidth',2)
title_string = strcat(['Meromorphic fct in the unit disk from boundary values, $f(z) = 100\sec($' num2str(beta2) '$z)\csc($' num2str(beta2) '$z)$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
xlabel('Re','Interpreter','LaTex')
ylabel('Im','Interpreter','LaTex')
set(gca,'FontSize',20);
plot(polesrat2,'r*','LineWidth',2)
leg{1} = ['Max Error of Ratdisk Approximant of type $($' num2str(mplot2-1) '$,$' num2str(mplot2-1) '$)$']';
leg{2} = 'Unit Sphere';
leg{3} = 'Poles of Ratdisk Approximant';
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
hold off