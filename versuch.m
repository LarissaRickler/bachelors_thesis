fun = @(x) sin(x);
Z = linspace(0,8,10000);

[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-13,10);


figure 
fplot(fun,[0,8],'-','Color',[1 0 0],'LineWidth',5)
hold on 
fplot(r,[0,8],'-','Color',[0 0 1],'LineWidth',3) % AAA-App

a = 1.000001; b = 0.000001;
fun2 = @(x) a*fun(x) + b;
[r2,pol2,res2,zer2,z2,f2,w2,errvec2] = aaa(fun2(Z),Z,1e-13,10);
r3 =@(x) a*r(x) + b;

figure 
fplot(fun2,[0,8],'-','Color',[1 0 0],'LineWidth',5)
hold on 
fplot(r2,[0,8],'-','Color',[0 0 1],'LineWidth',3) % AAA-App
fplot(r3,[0,8],'--','Color',[0 1 1],'LineWidth',3)
