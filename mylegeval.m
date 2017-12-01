function y = legpolyval(c, x, dom)
% Evaluate mapped-Legendre expansions on domain dom.
% (c = coefficients, x = evaluation points.)

% Scale values to [-1 1];
if ( nargin > 2 ), x = 2*x/dom(2) - 1; end

% Clenshaw scheme:
bk1 = 0*x; 
bk2 = bk1;
n = size(c,1)-1;
for k = n:-1:1
    bk = c(k+1) + (2*k+1)/(k+1)*x.*bk1 - (k+1)/(k+2)*bk2;
    bk2 = bk1;
    bk1 = bk;
end
y = c(1) + x.*bk1 - .5*bk2;

end