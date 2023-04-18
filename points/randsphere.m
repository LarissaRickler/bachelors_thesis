function Z=randsphere(mx,my,rc,n)
% INPUT: x1,y1 = coordinates of center
% rc = radius of circle
%
% OUTPUT: x,y = coordinates of random point on sphere
Z = zeros(n,1);
for k = 1:n
    a = 2*pi*rand;
    Z(k) = ((rc)*cos(a)+mx) + 1i*((rc)*sin(a)+my);
end
end