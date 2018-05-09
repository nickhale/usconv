% Example 1 from VIDECT paper
% Nick Hale - 2018

%% Set-up problem:

% Domain:
T = 1;
d = [0 T];

% Parameters:
a = 100;
b = sqrt(a^2-2*a+5)/2;

% Kernel:
K = @(t) exp(-t);

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

% % Faster / more accurate computation of Legendre coefficients:
% K_ = @(N) [ones(1, N) ; speye(N-1,N)*(D(N)+S12(N))]\[exp(-1); zeros(N-1,1)];
% V = @(N) Vmat(N, speye(N, 2*N)*K_(2*N), d, 1e-20);

tt = linspace(0, T, 1000);
err = [];
err2 = [];
NN = [2:100 980:1000];
y_ = NaN;
for N = [NN 60]
    fprintf('%d ', N);

    % The integro-differential operator:
    A = D(N) + S12(N)*(a*I(N) - V(N));
    
    % Boundary conditions:
    B = ones(1,N);
    B(1,2:2:end) = -1;      % Dirichlet Left
    gam = 1;
    
    % Almost-banded operator:
    A = [B ; A(1:N-1,:)];
    
    % RHS:
    f_ = zeros(N,1);
    rhs = [gam ; f_(1:N-1)];
    
    % Approximate coefficients of solution:
%     y_ = mysolve(A, rhs, 1); % Schur factorisation. 
    y_ = A\rhs;            % Backslash.
    
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
% print -depsc2 example1_spy

% Solution:
figure(2)
plot(tt, y(tt), 'LineWidth', 3);
ylim([0 1])
% print -depsc2 example1sol

% Error:
figure(3)
semilogy(NN, err(NN), '-', 'LineWidth', 3), shg
axis([0 NN(end) 1e-16, 1e1]), grid on
h = breakxaxis([100 980]);
set(h.leftAxes, 'XTick', [0 20 40 60 80 100])
set(h.rightAxes, 'XTick', 1000, 'fontsize', 9.5)
% print -depsc2 example1_err

% Solution (linear and log)
figure(4)
[ax, h1, h2] = plotyy(tt, y(tt), tt, (y(tt)), @plot, @semilogy);
set(h1, 'LineWidth', 3);
set(h2, 'LineWidth', 3, 'LineStyle', ':');
set(ax(1), 'ylim', [0 1])
grid on
legend('plot(t, y)', 'semilogy(t, y)')
% print -depsc2 example1_logsol

% Align the figures for display:
alignfigs % (http://github.com/nickhale/alignfigs)
