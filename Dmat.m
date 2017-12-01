function D = Dmat(N, m, dom)
% NxN Legendre-US differentiation matrix of order m on domain dom.

% Default to m = 1;
if ( nargin == 1 ), m = 1; end
% Scale to domain:
if ( nargin < 3 ), const = 1; else, const = 2/diff(dom); end
% Construct operator:
const = (2*const)^m*gamma(m+.5)/sqrt(pi);
D = const*spdiags(ones(N, 1), m, N, N);

end
