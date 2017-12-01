function V = Vmat(N, K, dom, tol)
% NxN Legendre-US Voterra operator matrix with kernel K on the domain dom.
%
% K may be a function handle or specified by it's Legendre coefficients on 
% [0 dom(2)]. In either casecoefficients less than tol in magnitude are set to
% zero.
%
% This code is somewhat optimised for spead and accuracy. A more readable
% version is provided below, but it is not recommended that this be used for
% computation (in particular, the readable version will be unstable for even
% moderate values of N.)

    if ( nargin < 3 )
        dom = [0 1];
    end
    if ( nargin < 4 )
        tol = 100*eps;
    end
    
    if ( isa(K, 'function_handle') )
        K = mylegcoeffs(K, N, dom);
    end
    
    % Ensure K has length N:
    K = [K ; zeros(N-length(K),1)];
    if ( length(K) > N )
        K = K(1:N);
    end
    
    % Zero small coefficients:
    Kmax = norm(K, inf);
    K(abs(K) < Kmax*tol) = 0;
    
    % Uncomment the line below to use (unstable but) readable version
%     V = Vmat_readable(N, K, dom); return 
        
    % S represents multiplication by 1/z in spherical Bessel space:
    e = [[1 ; 1./(2*(1:(N-1)).'+1)], [1 ; zeros(N-1, 1)], -1./(2*(0:N-1).'+1)];
    S = spdiags(e, -1:1, N, N);

    % Initialise scl:
    scl = 1./(2*(1:N).'-1); scl(2:2:end) = -scl(2:2:end); scl(1) = 0;
    
    % Initialise V:
    M = find(K~=0, 1, 'last');
    if ( M < N/2 ) 
        Vlower = sparse([], [], [], N, N, N*M);
        Vupper = Vlower;
        K = sparse(K);
    else
        Vlower = zeros(N);
        Vupper = Vlower;
    end
    
    % First column of B:
    v = diff(dom)/2*S*K; 
    vOld = v;
    Vlower(:,1) = v;
    % First row of B:
    Vupper(:,1) = v.*scl;

    % The scalar case is trivial:
    if ( N == 1 )
        V = full(Vlower + Vupper.');
        return
    end
    
    % Second column of B:
    v = S*v - v; v(1) = 0;
    Vlower(:,2) = v;
    % Second row of B:
    scl = -scl*3; scl(2) = 0;
    Vupper(:,2) = v.*scl;
    
    % Loop over remaining columns using recurrence:
    for n = 3:N
        % Recurrence
        vNew = (2*n-3) * (S * v) + vOld;  
        vOld = v; 
        v = vNew;
        v([n-2,n-1]) = 0;
        % Add to Bl
        Vlower(:,n) = v;
        % Recurrence is unstable for j < k. Correct for upper-tri part:
        scl = -scl*((n-.5)/(n-1.5)); scl(n) = 0;
        Vupper(:,n) = v.*scl;
    end
    
    % Recombine lower and upper for the result:
    V = Vlower + Vupper.';
    
end

function W = Vmat_readable(N, K, dom)
% A more human readable version of the above, but less efficient and less
% accurate. (In particular, the recurrence is unstable for j < k when n >> 1.

% Initialise (pad with an extra row and column for simplicity)
W = zeros(N+1,N+1); K(end+1) = 0;

% Useful values:
scl = diff(dom)/2;
j = (1:N-1)';
n = (0:N)';

% W_00 entry:
W(1,1)   = scl * ( K(1) - K(2)/3 );
% First column:
W(j+1,1) = scl * ( K(j)./(2*j-1) - K(j+2)./(2*j+3) );

% First row:
W(1,:)       =  W(:,1)./(2*n+1);
W(1,2:2:end) = -W(1,2:2:end);

% W_01 entry:
W(1,2)   = -W(2,1)/3;
% Second column:
W(j+1,2) =  W(j,1)./(2*j-1) - ...
            W(j+1,1) - ...
            W(j+2,1)./(2*j+3);

% Recurrence relation:
for n = 1:N-2
    W(j+1,n+2) = -(2*n+1)./(2*j+3).*W(j+2,n+1) + ...
                  (2*n+1)./(2*j-1).*W(j,n+1) + ...
                   W(j+1,n);
end

% Trim padding:
W = W(1:N, 1:N);

end


