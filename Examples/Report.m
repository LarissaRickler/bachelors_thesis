%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   Sawtooth wave  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sawtooth wave function
g = @(x) x - floor(x);

% interval
a = 0; b = 2;

% AAA-interpolant with fake nodes
nfn = 1000; % number of samplepoints
xfn = linspace(a,b,nfn); % equispaced sample points
k = 50; % shifting parameter
S = @(x) sGibbs((a + b)/2,2,k,x); % fake nodes map
G = g(xfn); % function values at samplepoints
[rfn,polfn,resfn,zerfn,zfn,ffn,wfn,errvecfn] = aaa_mod(G,xfn,S,1e-13,100); % interpolant

% (sawtooth wave) - (AAA interpolant with fake nodes)
errfunfn = @(x) g(x) - rfn(x); 
% alternativ error
xx = linspace(a,b,10000);
errintfn = g(xx) - rfn(xx);

% plot
fplot(g,[a,b],'-','Color',[1 0 0],'LineWidth',5)
hold on 
fplot(rfn,[a,b],'-','Color',[0 0 1],'LineWidth',3)

leg{1} = 'Sawtooth Wave';
leg{2} = 'AAA-interpolant with fake nodes';
title_string = strcat('Sawtooth Wave');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off
 
figure
fplot(errfunfn,[a,b],'-','Color',[1 0 0],'LineWidth',3)
hold on
leg{1} = 'error as function';
title_string2 = strcat('Error Computation of AAA-interpolant with fake nodes');
title(title_string2,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off 


figure
plot(xx,errintfn,'-','Color',[0 0 1],'LineWidth',3)
hold on
leg{1} = 'error of equidistant points';
title_string2 = strcat('Error Computation of AAA-interpolant with fake nodes');
title(title_string2,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off