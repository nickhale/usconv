% Example 2 from VIDECT paper

%% Set-up problem:

% Domain:
T = 1;
d = [0 T];

% Parameters:
a = 100;
b = sqrt(a^2-2*a+5)/2;

% Kernel:
K = @(t) exp(-t);

% RHS:
f = @(t) exp(-t)/b.* ...
    (exp(-(a-1)/2*T)*sinh(b*T) - exp(-(a-1)/2*t).*sinh(b*t));

% Exact solution
sol = @(t) exp(-(a+1)/2*t).*(cosh(b*t) - .5*(a-1)/b*sinh(b*t));

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

tt = linspace(0, T, 1000);
err = [];
err2 = [];
NN = [2:100 980:1000];
for N = [NN 60]
    fprintf('%d ', N);
    
    % The integro-differential operator:
    A = D2(N) + S32(N)*(a*D(N) + S12(N)*(F(N) - eye(N)));
    
    % Boundary conditions:
    B = ones(2,N);          % Dirichlet right
    B(1,2:2:end) = -1;      % Dirichlet Left
    gam = [1 ; sol(T)];
    
    % Almost-banded operator:
    A = [B ; A(1:N-2,:)];
    
    % RHS:
    f_ = S32(N)*S12(N)*mylegcoeffs(f, N, d);
    rhs = [gam ; f_(1:N-2)];
    
    % Approximate coefficients of solution:
    y_ = mysolve(A, rhs, 2); % Schur factorisation. 
%     y_ = A\rhs;            % Backslash.
    
    % Error:
    y = @(t) mylegeval(y_, t, d);
    err(N) = norm(y(tt) - sol(tt), inf);
    
end
fprintf('\n');

%%
close all

% Spy plot:
figure(1)
spy(A)
% print -depsc2 example2_spy

% Solution:
figure(2)
plot(tt, y(tt), 'LineWidth', 3);
% print -depsc2 example2sol

% Error:
figure(3)
semilogy(NN, err(NN), '-', 'LineWidth', 3), shg
axis([0 NN(end) 1e-16, 1e1]), grid on
h = breakxaxis([100 980]);
set(h.leftAxes, 'XTick', [0 20 40 60 80 100])
set(h.rightAxes, 'XTick', 1000)
% print -depsc2 example2_err

% Solution (linear and log)
figure(4)
[ax, h1, h2] = plotyy(tt, y(tt), tt, (y(tt)), @plot, @semilogy);
set(h1, 'LineWidth', 3);
set(h2, 'LineWidth', 3, 'LineStyle', ':');
set(ax(1), 'ylim', [0 1])
grid on
legend('plot(t, y)', 'semilogy(t, y)')
% print -depsc2 example2_logsol

% Align the figures for display:
alignfigs % (http://github.com/nickhale/alignfigs)









