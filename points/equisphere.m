function Z = equisphere(mx,my,rc,n)
% INPUT: m1,m1 = coordinates of center
% rc = radius of circle
% n = number of points
%
% OUTPUT: Z = coordinates of n equispaced points on complex sphere
Z = zeros(n,1);
for k = 1:n
    a = (2*pi*k)/n;
    Z(k) = ((rc)*cos(a)+mx) + 1i*((rc)*sin(a)+my);
end
end