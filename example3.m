% Example 3 from VIDECT paper

%% Set-up problem:

% Domain:
T = 1;
d = [0 T];

% Parameters:
a = .1;
b = .2;
r = sqrt(a^2+b^2);

% Kernel:
K = @(t) exp(-.5*t.^2/b^2);

% RHS:
f = @(t) a*b*sqrt(pi/2)/r*exp(-.5*(t/r).^2).* ...
    (erf(a/(b*r)*t/sqrt(2)) + erf((r^2-a^2*t)/(a*b*r*sqrt(2))));

% Exact solution
sol = @(t) exp(-.5/a^2*t.^2);

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
% Multiplication operator:
M = @(N) Mmat(N, @(t) t, d, 3/2); % t in C^{3/2}

err = [];
tt = linspace(0, T, 1000);
NN = [2:100 980:1000];
for N = [NN 60]
    fprintf('%d ', N);
    
    % The integro-differential operator:
    A = a^2*D2(N) + S32(N)*( M(N)*D(N) + S12(N)*(I(N) + F(N)) );
    
    % Boundary conditions:
    B = ones(2,N); 
    B(1,2:2:end) = -1;  % Dirichlet
    B(2,2:end) = 0;     % Mean value
    gam = [1 ; sqrt(pi/2)*a*erf(1/(sqrt(2)*a))];
    
    % Almost-banded operator:
    A = [B ; A(1:N-2,:)];
    
    % RHS:
    f_ = S32(N)*S12(N)*mylegcoeffs(f, N, d);
    rhs = [gam ; f_(1:N-2)];
    
    % Approximate coefficients of solution:
    y_ = mysolve(A, rhs, 2); % Schur factorisation. 
%     y_ = A\rhs;            % Backslash.
    
    % Approximate infinity norm err
    y = @(t) mylegeval(y_, t, d);
    err(N) = norm(y(tt) - sol(tt), inf);
    
end
fprintf('\n');

%% Plotting:

% Spy plot:
figure(1)
spy(A);
% print -depsc2 example3_spy

% Solution:
figure(2)
plot(tt, y(tt), 'LineWidth', 3);
% hold on, plot(tt, sol(tt), '--', 'LineWidth', 2), hold off
ylim([0 1])
% print -depsc2 example3_sol

% Error:
figure(3)
semilogy(NN, err(NN), '-', 'LineWidth', 3), shg
axis([0 NN(end) 1e-16, 1e1]); grid on
h = breakxaxis([100 980]);
set(h.leftAxes, 'XTick', [0 20 40 60 80 100]);
set(h.rightAxes, 'XTick', 1000);
% print -depsc2 example3_err

% Solution (linear and log)
figure(4)
[ax, h1, h2] = plotyy(tt, y(tt), tt, (y(tt)), @plot, @semilogy);
set(h1, 'LineWidth', 3);
set(h2, 'LineWidth', 3, 'LineStyle', ':');
set(ax(1), 'ylim', [0 1]);
set(ax(2), 'ylim', [1e-15 1e0]);
grid on
legend('plot(t, y)', 'semilogy(t, y)');
% print -depsc2 example3_logsol

% Align the figures for display:
alignfigs % (http://github.com/nickhale/alignfigs)








