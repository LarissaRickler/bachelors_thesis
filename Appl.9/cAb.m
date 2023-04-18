function f = cAb(c,A,b,x)
dim = length(c);
f = zeros(length(x),1);
for k = 1:length(x)
    f(k) = c*((x(k)*eye(dim)-A)^(-1))*b;
end
end