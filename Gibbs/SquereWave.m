%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Square Wave %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% interval [a, b] 
a = 0; b = 2; 
global L 
L = (b - a)/2;

% square wave function
f = @(x) 2.*(heaviside(x/L) - heaviside(x/L - 1)) - 1;

% fourier series of degree nf
nf = 100; % degree of fourier series
t = @(x) fourierSquareWave(x, nf); % interpolant

% AAA-Interpolant 
nr = 1000; % number of sample points
xr = linspace(a,b,nr); % equispaced sample points
% xr = a + (b - a).*rand(nr, 1); % uniformly distributed sample points
[r,pol,res,zer,z,fj,w,errvec] =  aaa(f,xr,1e-13,10); % interpolant

% AAA-Interpolant with fake nodes
nfn = 1000; % number of sample points
xrfn = linspace(a,b,nfn); % equispaced sample points
% xrfn = a + (b - a).*rand(nfn, 1); % uniformly distributed sample points
k = 500; % shifting parameter
S = @(x) sGibbs(L,2,k,x); % fake nodes map
F = f(xrfn); % function values at sample points 
[rfn,polfn,resfn,zerfn,zfn,fjfn,wfn,errvecfn] = aaa_mod(F,xrfn,S,1e-13,10); % interpolant

% plot
fplot(f, [a, b], '-', 'Color', [1,0.6,0], 'LineWidth', 3)
hold on
fplot(t, [a, b], '-', 'Color', [1,0,0], 'LineWidth', 1)
fplot(r, [a, b], '-', 'Color', [0,0,1], 'LineWidth', 3)
fplot(rfn, [a, b], '--', 'Color', [1,0,1], 'LineWidth', 3)

leg{1} = 'Square Wave';
leg{2} = 'Fourier Series';
leg{3} = 'rational Interpolant';
leg{4} = 'rational Interpolant with fake nodes';
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
