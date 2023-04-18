function S = sGibbs(xi,d,k,xx)
% Input: xi = position of jump discontinuties
% d = magnitdues of jump discontinuties
% k = shifting parameter > 0
% xx = argument
%
% Output: S = value of fake-nodes-mapping at xx
n = length(xx);
S = zeros(n,1);
for i = 1:n
    S(i) = sGibbsHelper(xi,d,k,xx(i));
end
end

function S = sGibbsHelper(xi,d,k,xx)
alpha = k.*d(:);
m = length(d);
l = 0;
while l < m & xx >= xi(l+1)
    l = l + 1; 
end
if l ~= 0
    A = sum(alpha(1:l));
    S = xx + A;
else
    S = xx;
end
end