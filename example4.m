% Example 4 from VIDECT paper
% Nick Hale - 2018

%% Set-up problem:

% Domain:
T = 1;
d = [0, T];

% Parameters:
w = 20;
n = 3;
m = 2;
gam = [.5*(n==1) ; (n==2)*w/4];

% Kernel:
J = @(n,x) besselj(n,x);
K = @(t) w*J(m,w*t);

% RHS:
f = @(t) J(m+n,w*t) + .5./t.^2.*((n-1)*(n-2)*J(n-1,w*t)+(n+1)*(n+2)*J(n+1,w*t));

% Exact solution
sol = @(t) n*J(n, w*t)./(w.*t);

%% New method:

% Differential operators:
D = @(N) Dmat(N, 1, d);
D2 = @(N) Dmat(N, 2, d);
% Integral operators:
V = @(N) Vmat(N, K, d);
F = @(N) Fmat(N, K, d);
% Conversion operators:
I = @(N) speye(N);
S12 = @(N) Smat(N, 1/2);
S32 = @(N) Smat(N, 3/2);

err = [];
tt = linspace(0, T, 1000); 
tt(1) = []; % Avoid removable singularity at t = 0.

NN = [2:100 990:1000];
for N = [NN,  60]
    fprintf('%d ', N);
    
    % The integro-differential operator:
    A = D2(N) + S32(N)*S12(N)*( w^2*I(N) + V(N) );
    
    % Boundary conditions:
    B = ones(2,N); 
    B(1,2:2:end) = -1;                           % Dirichlet left
    nn = 0:N-1; B(2,:) = T*nn.*(nn+1).*(-1).^nn; % Neumann left
    
    % RHS:
    f_ = S32(N)*S12(N)*mylegcoeffs(f, N, d);
    rhs = [gam ; f_(1:N-2)];     
    
    % Almost-banded operator:
    A = [B ; A(1:N-2,:)];
    
    % Approximate coefficients of solution:
%     y_ = mysolve(A, rhs, 2, 1); % Schur factorisation. 
    y_ = A\rhs;            % Backslash.

    % Approximate infinity norm error:
    y = @(t) mylegeval(y_, t, d);
    err(N) = norm(y(tt) - sol(tt), inf);
    
end
fprintf('\n');

%% Plotting:
close all
 
% Spy plot:
figure(1)
spy(A)
% print -depsc2 example4_spy

% Solution:
figure(2)
plot(tt, y(tt), 'LineWidth', 3);
% print -depsc2 example4_sol

% Error:
figure(3)
semilogy(NN, err(NN), '-', 'LineWidth', 3), shg
axis([0 NN(end) 1e-16, 1e1]), grid on
hold on
h = breakxaxis([50 990]);
set(h.leftAxes, 'XTick', [0 10 20 30 40 50]);
set(h.rightAxes, 'XTick', 1000, 'FontSize', 9.5);
print -depsc2 example4_err

% Align the figures for display:
alignfigs % (http://github.com/nickhale/alignfigs)
