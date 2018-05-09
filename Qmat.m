function Q = Qmat(N, dom)
% NxN Legendre-US indefinite integration matrix on domain dom.

% Scale to domain:
if ( nargin < 2 ), const = 1; else, const = diff(dom)/2; end

% Construct operator:
v = 1./(1:2:2*N-1)';
Q = const*spdiags([v, -v], [-1 1], N, N);

if ( N > 1 )
    Q(1) = Q(2);
end

end