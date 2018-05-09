% Example 4 from VIDECT paper

%% Set-up problem:

% Domain:
T = 1;
d = [0, T];

% Parameters:
w = 20;
n = 3;
m = 2;
Eps = .001;
gam = [.5*(n==1) ; (n==2)*w/4];

% % Kernel:
J = @(n,x) besselj(n,x);
K = @(t) w*J(m,w*t);
K2 = chebfun2(@(x,y) K(x-y), [0 T 0 T]);
g = chebcoeffs(K2.rows);
h = chebcoeffs(K2.cols);
R = size(g, 2);

% RHS:
f = @(t) J(m+n,w*t) + .5./t.^2.*((n-1)*(n-2)*J(n-1,w*t)+(n+1)*(n+2)*J(n+1,w*t));

% Exact solution
sol = [];

%% Approxfun method (Chebyshev basis):

% Differential operators:
D2 = @(N) ultraS.diffmat(N, 2)*(2/T)^2;
% Integral operators:
Q = @(N) Qmat_cheb(N, d);
% Conversion operators:
I = @(N) speye(N);
S = @(N) ultraS.convertmat(N,0,1);

err = [];
tt = linspace(d(1), d(2), 1000); 
tt(1) = []; % Avoid removable singularity at t = 0.
NN = [2000 2:400 900:1000];
NN = []; tic
for N = [NN,  400]
    fprintf('%d ', N);
    
    QN = Q(N);
    QN(abs(QN)<1e-10) = 0;
    V = sparse(0);
    for r = 1:R
        Mg = ultraS.multmat(N, g(:,r), 0);
        Mh = ultraS.multmat(N, h(:,r), 0);
        V = V + 1./K2.pivotValues(r)*Mg*QN*Mh;
    end
    
    % The integro-differential operator:
    A = Eps*D2(N) + S(N)*( w^2*I(N) + V );
    
    % Boundary conditions:
    B = ones(2,N); 
    B(1,2:2:end) = -1; % Dirichlet left
    nn = 0:N-1; 
    B(2,:) = nn.^2.*(-1).^nn*2/T; % Neumann left
    
    % RHS:
    f_ = S(N)*chebcoeffs(chebfun(f, d), N);
    rhs = [gam ; f_(1:N-2)];     
    
    % Almost-banded operator:
    A = [B ; A(1:N-2,:)];
    
    % Approximate coefficients of solution:
    y_ = A\rhs;            % Backslash.
    
    if ( isempty(sol) )
        % Use large N for 'true' solution.
        sol = chebfun(y_, d, 'coeffs');
    end

    % Approximate infinity norm error:
    y = chebfun(y_, d, 'coeffs');
    err(N) = norm(y(tt) - sol(tt), inf);
    
end
toc, return
fprintf('\n');

%% Plotting:

% Spy plot:
figure(5)
spy(A)
print -depsc2 example4b_approxfun_cheb_spy

% Error:
figure(3)
semilogy(NN, err(NN), '-', 'LineWidth', 3), shg
axis([0 NN(end) 1e-16, 1e1]), grid on

% Align the figures for display:
alignfigs % (http://github.com/nickhale/alignfigs)

function Q = Qmat_cheb(N, dom)
v = (2:2:2*N)';
v2 = -[0 ; 0 ; v(1:end-2)];
Q = spdiags(1./[v, v2], [-1 1], N, N);
v3 = (0:N-1).^2-1;
v3(1:2:end) = -v3(1:2:end);
Q(1,:) = 1./v3;
Q(1:2,1:2) = [1 -.25 ; 1 0];
Q = Q*diff(dom)/2;
end







