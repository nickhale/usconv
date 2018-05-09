function F = Fmat(N, K, varargin)
% NxN Legendre-US Fredholm operator matrix with kernel K on the domain dom.
%
% K may be a function handle, in which case the calling sequence should be 
% FMAT(N, K, DOM, TOL) (where coefficients less than tol in magnitude are set to
% zero).
%
% Alternatively K may be specified by its Legendre coefficients on [0 dom(2)], 
% in which case the calling sequence is FMAT(N, K, KHAT, DOM, TOL), where KHAT
% are the Legendre coefficients of the kernel on the inteval [-dom(2), 0].

% Parse the inputs:
if ( isnumeric(K) )
    Khat = varargin{1};
    varargin(1) = [];
else
    Khat = @(t) K(-t);
end

% Volterra opertor:
V1 = Vmat(N, K, varargin{:});

% 'flip' matrix:
e = ones(N,1); 
e(2:2:end) = -1; 
I1 = spdiags(e, 0, N, N);

% 'flipped' Volterra operator:
V2 = I1*Vmat(N, Khat, varargin{:})*I1;

% Fredholm operator:
F = V1 + V2;
    
end