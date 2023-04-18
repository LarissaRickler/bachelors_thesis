%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Condition Number Lambda %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% and Lambda C %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% 1 Analytic fct in the unit disk
% tan(x)
disp('1 Analytic fct in the unit disk')
M = 128; % # sample points, similar results for m >> 128
m = 7; % type (m-1,m-1) paper m <= 9
Z = equisphere(0,0,1,M); % M equispaced points in complex unit sphere 
fun = @(x) tan(x);
[r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-13,m);
% Grid for Lambda
G1 = disk_grid(100,1,[0,0],1);
G = G1(1,:) + 1i*G1(2,:);


% % 2 Analytic fct in unit disk with nearby brance points
% % log(1.1 - x)
% disp('2 Analytic fct with nearby brance points')
% M = 256; % # sample points, paper 256
% m = 12; % type(m-1,m-1), paper: AAA m <= 16, ratdisk m <= 20
% Z = equisphere(0,0,1,M); % M equispaced points in complex unit circle  
% fun = @(x) log(1.1 -x);
% [r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-13,m);
% % Grid for Lambda
% G1 = disk_grid(100,1,[0,0],1);
% G = G1(1,:) + 1i*G1(2,:);


% % 3 Meromorphic fct in the unit disk from boundary values
% % tan(beta * x) with beta = 256 in the unit disk
% disp('3 Meromorphic fct in the unit disk from boundary values')
% M = 1000; % # sample points, paper 1000
% m = 12; % type(m-1,m-1), paper: beta/rage(m): 4/(0,14), 16/(0,28), 64/(0,49), 256/(0,62)
% Z = equisphere(0,0,1,M); % M points on unit sphere
% beta = 4; % beta = 4,16,64,256
% fun = @(x) tan(beta*x);
% [r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-13,m);
% % Grid for Lambda
% G1 = disk_grid(100,1,[0,0],1);
% G = G1(1,:) + 1i*G1(2,:);


% % 4 Meromorphic fct in the unit disk from bound and int values
% % tan(beta * x) with beta = 256 in the unit disk from boundary and
% % interior values
% disp('4 Meromorphic fct in the unit disk from bound and int values')
% M1 = 1000; % # sample points on unit sphere, paper 1000
% M2 = 3000; % # sample points in unit disk, paper 3000
% m = 12; % type(m-1,m-1), paper: beta/rage(m): circa 4/(0,14), 16/(0,28), 64/(0,65), 256/(0,200)
% Z1 = equisphere(0,0,1,M1); % M1 equispaced points in unit sphere
% Z2 = randdisk(0,0,1,M2); % M2 random points in unit disk
% Z = [Z1; Z2];
% beta = 4; % beta = 4,16,64,256 
% fun = @(x) tan(beta*x);
% [r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-13,m);
% % Grid for Lambda
% G1 = disk_grid(100,1,[0,0],1);
% G = G1(1,:) + 1i*G1(2,:);


% % 5 Appr. in connected domains
% % 1/bessel
% disp('5 Appr. in connected domains')
% M = 2000; % # sample points, paper 2000
% m = 13; %max type (m-1,m-1), paper m = 13
% fun = @(x) 1./besselj(0,x); %1/J_0 with J_0 bessel fct
% Z = 10*rand(M,1) + 1i*(-1+2*rand(M,1)); % M random points in complex plane (0,10)+i(-1,1)
% [r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-13,m);
% % Grid for Lambda
% [GX,GY] = meshgrid(linspace(0,10,1000),linspace(-1,1,200));
% G1 = GX + 1i*GY;
% G = G1(:)';


% % 6 Appr. in disconnected domains
% % sign(Re(x))
% disp('6 Appr. in disconnected domains')
% m = 51; % type (m-1,m-1), paper 45-> actually type(43,43), convergence at 51
% M = 1000; % # points in each section, paper 1000
% Z1 = equisphere(1.5,0,1,M); % M equispaced on unit circle
% Z21 = linspace(-2.5,-0.5,M/4 + 1);
% Z21 = Z21(1:end-1);
% Z22 = linspace(-0.5,-2.5,M/4 + 1);
% Z22 = Z22(1:end-1);
% Z23 = linspace(1,-1,251);
% Z23 = Z23(1:end-1);
% Z24 = linspace(-1,1,251);
% Z24 = Z24(1:end-1);
% Z2 = [ Z21 + 1i, Z22 - 1i, -0.5 + 1i*Z23, -2.5 + 1i*Z24]; % M equispaced around square
% Z = [Z1;Z2'];
% % % M1 sample points on the unit sphere (equispaced)
% % M1 = 1000; % # sample points in sphere section
% % Z2 = equisphere(1.5,0,1,n1);
% % % M2^2 sample points in rectangle (equally spaced)
% % M2 = 32;
% % [X,Y] = meshgrid(linspace(-2.5,-0.5,M2),linspace(-1,1,M2));
% % Z1 = X + 1i*Y;
% % Z1 = Z1(:);
% % Z = [Z1;Z2];
% fun = @(x) sign(real(x));
% [r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-13,m);
% % Grid for Lambda
% [GX,GY] = meshgrid(linspace(-2.5,-0.5,173),linspace(-1,1,173));
% G1 = GX + 1i*GY;
% G2 = G1(:)';
% G3 = disk_grid(100,1,[1.5,0],10);
% G4 = (G3(1,:) + 1i*G3(2,:));
% G = [G2, G4];


% % 9 Clamped beam model form Chahlaoui and Van Dooren
% % fun(z) = c'(zI-A)^(-1)b (a rat. fun. of type(348,348) 
% disp('9 Clamped beam model form Chahlaoui and Van Dooren')
% load('matlab.mat')
% M = 500; % 2*n sample points
% m = 3; % type (m-1,m-1), paper between 0 and 46
% Z1 = 1i*logspace(-2,2,M); % M logarithmically spaced points from 10?2 i to 102 i
% Z2 = -Z1; % complex conjugates
% dim = 348;
% Z = [Z1 Z2];
% fun = @(x) cAb(C,A,B,x);
% [r,pol,res,zer,z,f,w,errvec] = aaa(fun(Z),Z,1e-13,m);
% % Grid for Lambda
% G = 1i*linspace(-100,100,20000);



numL = @(x) sum(abs(w)./abs(x-z)); % numerator of Lebesgue fct
denomL = @(x) abs(sum(w./(x-z))); % denominator of Lebesgue fct
Lambdaf = numL(G)./denomL(G); % Lebesgue fct on points G
Lambda = max(Lambdaf); % "Lebesgue constant" 

C = [];
Zm = setdiff(Z,z);
for k = 1:length(z)
C = [C 1./(Zm-z(k))];
end
LambdaC = cond(C, 2);




function f = cAb(c,A,b,x)
f = [];
for k = 1:length(x)
    f = [f; c*((x(k)*eye(348)-A)^(-1))*b];
end
end