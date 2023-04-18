%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Condition Number Lambda %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% and Lambda C %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all


% % 7 Approximation of |x| on [-1,1]
% % way one
% disp('7 Appr. of |x| on [-1,1], way one')
% M = 200000; % # samplepoints in [-1,1]
% m = 7; % max type (m-1,m-1), paper between 0 and 26
% fun = @(x) abs(x);
% Z = linspace(-1,1,M); % M equi spaced points in [-1,1]
% [r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-13,m);
% % Grid for Lambda
% G = Z;


% % 7 Approximation of |x| on [-1,1]
% % way two: approximate sqrt(x) on [0,1]
% disp('7 Appr. of |x| on [-1,1], way two')
% M = 100000; % # samplepoints in [-1,1]
% m = 80; % max type (m-1,m-1), (circa half of 'way one' for equal good approx.)
% fun1 = @(x) sqrt(x);
% Z = linspace(0,1,M); % M equi spaced points in [0,1]
% [r1,pol,res,zer,z,f,w,errvec] = aaa(fun1(Z),Z,1e-13,m);
% % Grid for Lambda
% G = [-Z,Z];
% % convert
% fun = @(x) fun1(x.^2);
% r = @(x) r1(x.^2);


% 8 Appr. of exp(x) in (-inf,0]
% Nicht sinnvoll, evtl wurde im Paper möbius trafo genutzt
disp('8 Appr. of exp(x) in (-inf,0]\n')
M = 4000; % # samplepoints
m = 7; % max type (m-1,m-1), paper between 0 and 14 
fun = @(x) exp(x);
Z = -logspace(-4,3,M); % M points logarithmically spaced between 10^-4 and 10^3
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-12,m);
% Grid for Lambda
G = -linspace(10^4,-10^(-3),100000);


numL = @(x) sum(abs(w)./abs(x-z)); % numerator of Lebesgue fct
denomL = @(x) abs(sum(w./(x-z))); % denominator of Lebesgue fct
Lambdaf = numL(G)./denomL(G); % Lebesgue fct on points G
Lambda = max(Lambdaf); % "Lebesgue constant" 

C = [];
Zm = setdiff(Z,z);
for k = 1:length(z)
C = [C; 1./(Zm-z(k))];
end
LambdaC = cond(C, 2);