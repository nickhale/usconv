function M = Mmat(N, f, dom, lam, tol)
% Construct the NxN US multiplication matrix for the function f on the domain
% dom in the space C^(lam). f may be a function handle or vector of coefficients
% in C^(lam). Coefficients below tol (default = 100*eps) are set to zero. If not
% given then dom defaults to [-1 1] and lam to 0.5 (Legendre).

% Default to [-1 1]:
if ( nargin < 3 )
    dom = [-1 1]; % Only needed if f is given as a function handle.
end
% Default to Legendre basis:
if ( nargin < 4 )
    lam = .5;
end
% Default tolerance:
if ( nargin < 5 )
    tol = 100*eps;
end

% C^(lam) coeffs:
if ( isa(f, 'function_handle') )
    c = mylegcoeffs(f, N, dom);
    c = ultra2ultra(c, .5, lam);                 
else
    c = f(:);
    if ( length(c) > N )
        c = c(1:N);
    end
end

% Zero entries in c below tolerance:
c(abs(c) < tol) = 0;

% Trivial cases:
if ( isempty(c) ), M = 0; return, end
if ( ~any(c) ), M = speye(N); return, end

% Trim trailing zeros:
c = c(1:find(c, 1, 'last'));

% Construct Jacobi matrix:
nn = (0:N)';
J = spdiags(.5*[(nn+1)./(nn+lam), (nn+2*lam-1)./(nn+lam)], [-1, 1], N, N);

% Initialise recurrence:
Cm1 = sparse(0); 
C = speye(N);
M = c(1)*C;

% Recurrence relation:
for n = 0:length(c)-2
    Cp1 = (2*(n+lam)/(n+1))*(J*C) - ((n+2*lam-1)/(n+1))*Cm1;
    Cm1 = C; 
    C = Cp1;
    if ( c(n+2) ~= 0 )
        M = M + c(n+2)*C;
    end
end

end