% Convergence of application 6
% 6 Appr. in disconnected domains
% Square wave

clear all


fun = @(x) sign(sin((2.*pi.*real(x))./2));
tol = 1e-13; % tolerance 1e-13/1e-14


M = 1000; % Sample points per domain
Z = [];
for k = 1:5 % # domains
    Z = [Z; equisphere(k-0.5,0,0.2,M)];
end

G = Z;



AllmaxErrorAAA = [];


for m = 1:80%max type (m-1,m-1)
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,tol,m); 

errfunAAA = @(x) abs(fun(x)-r(x));

AllmaxErrorAAA = [AllmaxErrorAAA, max(errfunAAA(G))];
end


figure 
plot(AllmaxErrorAAA,'o-','Color',[0 0 1],'LineWidth',2)
hold on
leg{1} = 'AAA Approximant of type(m-1,m-1)';
title_string = strcat(['Appr. in disconnected domains, $f(z) = $sign$(\sin(\frac{2\pi}{2}$Re$(z)))$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
set(gca,'FontSize',20);
xlabel('$m$','Interpreter','LaTex')
ylabel('max error','Interpreter','LaTex')
set(gca, 'YScale', 'log')
grid on
hold off

% levels of error
ReAx = linspace(0,5,2000);
ImAx = linspace(-2,2,2000);
[RE,IM] = meshgrid(ReAx,ImAx);
errorCont = errfunAAA(RE + 1i*IM);



figure
contour(RE,IM,errorCont,logspace(-12,-1,12)','LineWidth',2);
set(gca,'ColorScale','log')
colormap(parula)
hcb = colorbar;
set(gca, 'clim', [10e-12 10e-1])
hold on 
title_string = strcat(['Appr. in disconnected domains, $f(z) = $sign$(\sin(\frac{2\pi}{2}$Re$(z)))$']);
title(title_string,'Interpreter','LaTex','FontSize',20);
xlabel('Re','Interpreter','LaTex')
ylabel('Im','Interpreter','LaTex')
set(gca,'FontSize',20);
plot(Z,'.','LineWidth',2,'Color','k')
plot(pol,'r*','LineWidth',2)
plot(zer,'g*','LineWidth',2)
leg{1} = 'Contours of Error';
leg{2} = 'Sample points';
leg{3} = 'Poles of AAA Approximant';
leg{4} = 'Zeros of AAA Approximant';
legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
grid on 
hold off
