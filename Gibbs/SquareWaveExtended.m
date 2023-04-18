%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Square Wave %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extended  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interval [a, b] 
a = 0; b = 2; 
global L 
L = (b - a)/2;

% square wave function
f = @(x) 2.*(heaviside(x/L) - heaviside(x/L - 1)) - 1;

% fourier series of degree nf
nf = 100; % degree of fourier series 
t = @(x) fourierSquareWave(x, nf); % fourier series of degree nf


% AAA-Interpolant 
nr = 1000; % numper of sample points
xr = linspace(a, b, nr); % equispaces sample points
[r,pol,res,zer,z,fj,w,errvec] =  aaa(f, xr,1e-13,100); % interpolant

% AAA-Interpolant with fake nodes
nfn = 1000; % numper of sample points
xrand = linspace(a,b,nfn); % equispaces sample points
k = 500; % shifting parameter
S = @(x) sGibbs(L, 2, k, x); % fake nodes map
F = f(xrand); % function values of sample points
[rfn,polfn,resfn,zerfn,zfn,fjfn,wfn,errvecfn] = aaa_mod(F,xrand,S,1e-13,100); % interpolant

% polynomialinterpolant of deg np with equispaced nodes
np = 10; % degree of interpolant
xp = linspace(a, b, np); % equispaced nodes
fp = f(xp); % function values at nodes
pe = @(x) polyval(polyfit(xp, fp, np), x); % interpolant

% polynomialinterpolant of deg nt with tschebyscheff nodes 
nt = 10; % degree of interpolant
xt = tschebyscheff_punkte(a, b, nt); % tschebyscheff nodes
ft = f(xt); % function values at nodes
pt = @(x) polyval(polyfit(xt, ft, nt), x); % interpolant

% plot
fplot(f, [a, b], '-', 'Color', [1,0.6,0], 'LineWidth', 3)
hold on
fplot(t, [a, b], '-', 'Color', [1,0,0], 'LineWidth', 1)
fplot(r, [a, b], '-', 'Color', [0,0,1], 'LineWidth', 1)
fplot(pe, [a, b], '-', 'Color', [0,1,0], 'LineWidth', 1)
fplot(pt, [a, b], '-', 'Color', [1,1,0], 'LineWidth', 1)
fplot(rfn, [a, b], '-', 'Color', [0,0,0], 'LineWidth', 1)

leg{1} = 'Square wave';
leg{2} = 'Fourier serie';
leg{3} = 'AAA-interpolant with equispaced sample points';
leg{4} = 'polynominalinterpolation with equispaced nodes';
leg{5} = 'polynominalinterpolation with tschebyscheff nodes';
leg{6} = 'AAA-interpolant with fake nodes mapped sample points';
title_string = strcat('Square Wave');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
hold off


function t = fourierSquareWave(x, n)
global L
t = 0;
for k = 1:2:n
    t = t + (1/k)*sin((k*pi.*x)./L);
end
t = t.*4./pi;
end
