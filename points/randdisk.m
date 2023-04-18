function Z=randdisk(mx,my,rc,n)
% INPUT: x1,y1 = coordinates of center
% rc = radius of circle
%
% OUTPUT: x,y = coordinates of random point in disk
Z = zeros(n,1);
for k = 1:n
    a = 2*pi*rand;
    r = sqrt(rand);
    Z(k) = ((rc*r)*cos(a)+mx) + 1i*((rc*r)*sin(a)+my);
end
end