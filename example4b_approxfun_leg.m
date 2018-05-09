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
pic
% % Kernel:
J = @(n,x) besselj(n,x);
K = @(t) w*J(m,w*t);
K2 = chebfun2(@(x,y) K(x-y), [0 T 0 T]);
g = legcoeffs(K2.rows);
h = legcoeffs(K2.cols);
R = size(g, 2);

% RHS:
f = @(t) J(m+n,w*t) + .5./t.^2.*((n-1)*(n-2)*J(n-1,w*t)+(n+1)*(n+2)*J(n+1,w*t));

% Exact solution
sol = [];

%% Approxfun method (Legendre Basis):

% Differential operators:
D2 = @(N) Dmat(N, 2, d);
% Integral operators:
V = @(N) Vmat(N, K, d);
F = @(N) Fmat(N, K, d);
Q = @(N) Qmat(N, d);
% Conversion operators:
I = @(N) speye(N);
S12 = @(N) Smat(N, 1/2);
S32 = @(N) Smat(N, 3/2);

err = [];
tt = linspace(d(1), d(2), 1000); 
tt(1) = []; % Avoid removable singularity at t = 0.
NN = [2000 2:400 900:1000];
NN = []; 
for N = [NN,  400] 
    fprintf('%d ', N);
    
    QN = Q(N);
%     V = sparse(0);
%     for r = 1:R
%         Mg = Mmat(N, g(:,r), d, .5);
%         Mh = Mmat(N, h(:,r), d, .5);
%         V = V + 1./K2.pivotValues(r)*Mg*QN*Mh;
%     end
    V = makeV(N, g, h, QN, 1./K2.pivotValues);
    
    
    % The integro-differential operator:
    A = Eps*D2(N) + S32(N)*S12(N)*( w^2*I(N) + V );
    
    % Boundary conditions:
    B = ones(2,N); 
    B(1,2:2:end) = -1; % Dirichlet left
    nn = 0:N-1; B(2,:) = T*nn.*(nn+1).*(-1).^nn; % Neumann left
    
    % RHS:
    f_ = S32(N)*S12(N)*mylegcoeffs(f, N, d);
    rhs = [gam ; f_(1:N-2)];     
    
    % Almost-banded operator:
    A = [B ; A(1:N-2,:)];
    
    % Approximate coefficients of solution:
%     y_ = mysolve(A, rhs, 2, 1); % Schur factorisation. 
    y_ = A\rhs;            % Backslash.
    
    if ( isempty(sol) )
        % Use large N for 'true' solution.
        sol = @(t) mylegeval(y_, t, d);
    end

    % Approximate infinity norm error:
    y = @(t) mylegeval(y_, t, d);
    err(N) = norm(y(tt) - sol(tt), inf);
    
end
poc, return
fprintf('\n');

%% Plotting:

% Spy plot:
figure(6)
spy(A)
% print -depsc2 example4b_approxfun_leg_spy

% Error:
figure(3)
semilogy(NN, err(NN), '-', 'LineWidth', 3), shg
axis([0 NN(end) 1e-16, 1e1]), grid on
h = breakxaxis([400 900]);
set(h.leftAxes, 'XTick', [0 100 200 300 400]);
set(h.rightAxes, 'XTick', 1000);

% Align the figures for display:
alignfigs % (http://github.com/nickhale/alignfigs)

%%

function V = makeV(N, g, h, Q, c)

% Zero entries in c below tolerance:
lam = .5;
R = length(c);


% Construct Jacobi matrix:
nn = (0:N)';
J = spdiags(.5*[(nn+1)./(nn+lam), (nn+2*lam-1)./(nn+lam)], [-1, 1], N, N);

% Initialise recurrence:
Cm1 = sparse(0); 
C = speye(N);
Mg = {};
Mh = {};
for r = 1:R
    Mg{r} = g(1,r)*C;
    Mh{r} = h(1,r)*C;
end

% Recurrence relation:
for n = 0:length(g)-2
    Cp1 = (2*(n+lam)/(n+1))*(J*C) - ((n+2*lam-1)/(n+1))*Cm1;
    Cm1 = C; 
    C = Cp1;
    for r = 1:R
        Mg{r} = Mg{r} + g(n+2,r)*C;
        Mh{r} = Mh{r} + h(n+2,r)*C;
    end
end

V = 0;
for r = 1:R
    V = V + c(r)*Mg{r}*Q*Mh{r};
end

end