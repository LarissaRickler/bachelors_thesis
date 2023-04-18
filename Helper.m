
% Sample points
figure % Sample points
plot(ZG,'o')
hold on 
title_string = strcat(['Sample points for  m = ' num2str(m)]);
title(title_string,'Interpreter','LaTex','FontSize',20);
hold off

% Support points AAA
figure % Support points AAA
plot(zG,'o')
hold on 
title_string = strcat(['Support points for  m = ' num2str(m)]);
title(title_string,'Interpreter','LaTex','FontSize',20);
hold off


% % poles inside the Disk
% intpol = [intpol, sum(abs(polesrat2)<1)];


% % Unit circle
% range = 0:pi/50:2*pi;
% xunit = cos(range);
% yunit = sin(range);


% % Poles and zeros AAA
% figure 
% plot(pol(1:end),'s')
% hold on 
% plot(zer(1:end),'d')
% title_string = strcat(['poles and zeros for  m = ' num2str(m)]);
% title(title_string,'Interpreter','LaTex','FontSize',20);
% leg{1} = 'Poles';
% leg{2} = 'Zeros';
% legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
% hold off


% % Poles Ratdisk
% polesrat1 = complex(polesrat1)
% figure 
% plot(polesrat1,'o')
% hold on 
% plot(xunit, yunit,'--','Color', [0 0 0],'LineWidth',1)
% title_string = strcat(['poles and zeros for  m = ' num2str(m)]);
% title(title_string,'Interpreter','LaTex','FontSize',20);
% hold off


% % Real part of AAA Approximant (no need) 
% rratreal = @(x) real(rrat(x)); 


% % Real part of AAA Approximant
% rreal = @(x) real(r(x)); 

% figure 
% fplot(fun,[0,5],'-','Color',[1 0 0],'LineWidth',5) % Function
% hold on 
% fplot(rreal,[0,5],'-','Color',[0 0 1],'LineWidth',3) % AAA-App
% % fplot(rratreal,[-1,1],'--','Color',[0 1 0],'LineWidth',3) % Ratdisk
% leg{1} = 'Function';
% leg{2} = 'AAA Approximant';
% % leg{3} = 'Ratdisk Approximant';
% legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
% set(gca,'FontSize',20);
% hold off



% % Calculation for best approximation for f(z) = log(z-1.1)
% c = 1;
% kappa = (c - sqrt(c^2 -1))^2;
% K = @(x) ellipke(x);
% R = exp( (pi * K( (1- kappa^2) ) ) / (4 * K(kappa^2)));
% R
% R^2



% display('min, max, mean of Lambda')
% min(LambdaAll)
% max(LambdaAll)
% mean(LambdaAll)
% figure
% plot(LambdaAll,'o')
% 
% display('min, max, mean of LambdaC')
% min(LambdaCAll)
% max(LambdaCAll)
% mean(LambdaCAll)
% figure
% plot(LambdaCAll,'o')
% 
% display('min, max, mean of LambdaCVF')
% min(LambdaCAllVF)
% max(LambdaCAllVF)
% mean(LambdaCAllVF)
% figure
% plot(LambdaCAllVF,'o')


% figure 
% plot(imag(G),abs(funG),'-','Color',[1 0 0],'LineWidth',5)
% hold on 
% plot(imag(G),abs(rG),'-','Color',[0 0 1],'LineWidth',3) % AAA-App
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% 
% 
% leg{1} = 'Function';
% leg{2} = 'AAA Approximant';
% title_string = strcat('Clamped beam model form Chahlaoui and Van Dooren');
% title(title_string,'Interpreter','LaTex','FontSize',20);
% legend(leg,'Interpreter','LaTex','FontSize',20,'Location','EastOutside');
% set(gca,'FontSize',20);
% hold off
