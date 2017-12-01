function c = mylegcoeffs(f, N, dom, tol)
% Approximate N mapped-Legendre coeffs of the function f on the domain dom.
% |coefficients| < tol are set to zero. If not given, tol = eps.

% Note we simply evaluate on a 5*N Chebyshev grid to get the coefficients.
% This is not foolproof.
   
if ( N == 0 )
    % N = 0 is trivial case.
    c = [];
    return
end

% Default to machine eps:
if ( nargin < 4 ), tol = 100*eps; end

% Evaluate on a larger Chebyshev grid:
x = chebpts(5*N, dom);
c = chebvals2legcoeffs(f(x), 1);
c = c(1:N);
   
% Set small coefficients to zero:
c(abs(c) < max(abs(c))*tol) = 0;

end

function x = chebpts(N, dom)
% First-kind Chebsyhev points on [-1 1]:
x = sin(pi*((-N+1:2:N-1)/(2*N))).';
% Scale to dom = [dom(1), dom(2)]
x = dom(2)*(x + 1)/2 + dom(1)*(1 - x)/2;
end