function [x] = tschebyscheff_punkte(a,b,n)

    % Tschebyscheff-Punkte auf [-1,1]
    x = cos((2*(0:n)+1)/(2*n+2)*pi);
    
    % Transformation auf [a,b]
    x = (b-a)/2 * x + (a+b)/2;

end