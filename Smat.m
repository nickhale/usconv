function S = Smat(N, lam)
% NxN US conversion matrix. S(N,lam): % C^{(lam)}(x) --> C^{(lam+1)}(x)
if ( lam == 0 )
    % Lam = case is special (Chebyshev T are not US polynomials)
    e = ones(N,1);
    S = spdiags(.5*[e,-e], [0,2], N, N);
    S(1) = 1;
else
    e = lam./((0:N-1)'+lam);
    S = spdiags([e,-e], [0,2], N, N);
end
end