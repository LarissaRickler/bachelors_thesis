% Convergence of application 2
% 2 Analytic fct with nearby brance points
% atan((1/1.05).*x)

clear all

fun = @(x) atan((1/1.05).*x);
M = 301; % # sample points, 301
tol = 1e-14; % 1e-14;  % 1e-13/1e-14

% Sample points
Z = equisphere(0,0,1,M); % M equispaced points in complex unit circle  
% Z = randdisk(0,0,1,M); % M random points in complex unit disk
% Z1 = disk_grid(9,1,[0,0],1); % equispaced points in complex unit disk
% Z = Z1(1,:) + 1i*Z1(2,:);    % "
% Z = Z';                      % "

% Grid 
G1 = disk_grid(100,1,[0,0],1); % equispaced points in complex unit disk
G = G1(1,:) + 1i*G1(2,:);      % "


AllmaxErrorAAA = [];
AllmaxErrorRat = [];
AllmaxErrorFit = [];

for m=1:40
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,tol,m); 


% Ratdisk
[rrat,arat,brat,murat,nurat,polesrat,residuesrat] = ratdisk(fun,m-1,m-1,M,tol);


errfunAAA = @(x) abs(fun(x)-r(x));
errfunrat  = @(x) abs(fun(x)-rrat(x));


AllmaxErrorAAA = [AllmaxErrorAAA, max(errfunAAA(G))];
AllmaxErrorRat = [AllmaxErrorRat, max(errfunrat(G))];

end


np = [0:40];
perror =  1.05.^(-2*np); % Error of polynomial interpolation


figure 
plot(AllmaxErrorAAA,'o-','Color',[0 0 1],'LineWidth',2)
hold on
plot(AllmaxErrorRat,'o-','Color',[0 1 0],'LineWidth',2)
plot(perror,'o-','Color',[1 0 1],'LineWidth',2)
leg{1} = 'AAA Approximant of type(m-1,m-1)';
leg{2} = 'Ratdisk Approximant of type(m-1,m-1)';
leg{3} = 'Polynomial of degree 2(m-1)';
title_string = strcat('Analytic fct with nearby brance points, $f(z) = \arctan(\frac{1}{1.05}z)$');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','LaTex')
ylabel('max error','Interpreter','LaTex')
set(gca, 'YScale', 'log')
grid on
hold off

figure
plot(pol,'s')
hold on
plot(zer(2:end),'d')
leg{1} = 'Poles of AAA Approximant of type $(30,30)$';
leg{2} = 'Zeros of AAA Approximant of type $(30,30)$';
title_string = strcat('Analytic fct with nearby brance points, $f(z) = \arctan(\frac{1}{1.05}z)$');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('Re','Interpreter','LaTex')
ylabel('Im','Interpreter','LaTex')
grid on
hold off

figure
plot(polesrat,'s')
hold on
leg{1} = 'Poles of Ratdisk Approximant of type $(39,39)$';
title_string = strcat('Analytic fct with nearby brance points, $f(z) = \arctan(\frac{1}{1.05}z)$');
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',15,'Location','EastOutside');
xlabel('Re','Interpreter','LaTex')
ylabel('Im','Interpreter','LaTex')
grid on
hold off
