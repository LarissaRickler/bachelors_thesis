%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   Sawtooth wave  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sawtooth wave function
g = @(x) x - floor(x);

% interval
a = 0; b = 2;

% AAA-interpolant
n = 1000; % number of samplepoints
x = linspace(a,b,n); % equispaced sample points
[r,pol,res,zer,z,f,w,errvec] = aaa(g,x,1e-13,100); % interpolant

% AAA-interpolant with fake nodes
nfn = 1000; % number of samplepoints
xfn = linspace(a,b,nfn); % equispaced sample points
k = 50; % shifting parameter
S = @(x) sGibbs((a + b)/2,2,k,x); % fake nodes map
G = g(xfn); % function values at samplepoints
[rfn,polfn,resfn,zerfn,zfn,ffn,wfn,errvecfn] = aaa_mod(G,xfn,S,1e-13,100); % interpolant

% error
errfun = @(x) g(x) - r(x); % does not work, maybe because of epsilon_mach
errfunfn = @(x) g(x) - rfn(x);
xx = linspace(a,b,10000);
errintnormal = g(xx) - r(xx);
errintfn = g(xx) - rfn(xx);

% plot
fplot(g,[a,b],'-','Color',[1 0 0],'LineWidth',5) % plot function
hold on 
fplot(r,[a,b],'-','Color',[0 1 0],'LineWidth',3) % plot AAA approximant
fplot(rfn,[a,b],'-','Color',[0 0 1],'LineWidth',3) % plot AAA approximant with fake nodes

leg{1} = 'Sawtooth Wave';
leg{2} = 'AAA-interpolant';
leg{3} = 'AAA-interpolant with fake nodes';
title_string = strcat('Sawtooth Wave');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off
 
figure % plot error
fplot(errfun,[a,b],'-','Color',[0 0 1],'LineWidth',3) 
% plot error fct AAA approximant
hold on 
fplot(errfunfn,[a,b],'-','Color',[1 0 0],'LineWidth',3)
% plot error fct AAA approximant with fake nodes

% figure
% plot(xx,errintnormal,'-','Color',[0 1 0],'LineWidth',3) 
% % plot error points AAA 
% hold on
% plot(xx,errintfn,'-','Color',[0 0 1],'LineWidth',3) 
% % plot error points AAA with fake nodes

leg{1} = 'error of AAA-interpolant';
leg{2} = 'error of AAA-interpolant with fake nodes';
title_string2 = strcat('Error Computation');
title(title_string2,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off