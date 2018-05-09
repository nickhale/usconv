% Example 3 from VIDECT paper using Driscoll's integral reformulation
% Nick Hale - 2018

%% Set-up problem:

% Domain:
T = 1;
d = [0 T];

% Parameters:
a = .1;
b = .2;
r = sqrt(a^2+b^2);
gam = [1 ; sqrt(pi/2)*a*erf(1/(sqrt(2)*a))];

% Kernel:
K = @(t) exp(-.5*t.^2/b^2);

% RHS:
f = @(t) a*b*sqrt(pi/2)/r*exp(-.5*(t/r).^2).* ...
    (erf(a/(b*r)*t/sqrt(2)) + erf((r^2-a^2*t)/(a*b*r*sqrt(2))));

% Exact solution
sol = @(t) exp(-.5/a^2*t.^2);
         
%% Discretise:

err = [];
tt = linspace(0, T, 1000);
NN = [2:100 980:1000];
for N = [NN 60]
    fprintf('%d ', N);

    % Operators:
    [t, w] = chebpts(N, d);
    I = eye(N);
    Q = cumsummat(N, d);
    T = diag(t);
    S = diag(w);
    [X, Y] = ndgrid(t);
    KK = K(X-Y);
    
    % Main operator:
    A = [a^2*I + T*Q + (I+KK*S)*Q^2, 2*t + KK*S*t ; 
         w*Q^2, w*t]; % <-- Boundary constraints
    
    rhs = [f(t) - 1 - KK*S*ones(N,1) ;  gam(2) - w*ones(N,1)];
    
    % Solve for v
    v = A\rhs;
    c1 = v(end);
    v = v(1:end-1);
    
    % Reconstruct solution y (chebyshev coefficients):
    y_ = Q*(Q*v) + c1*t + 1;
    
    % Approximate infinity norm err
    y = chebfun(y_, d);
    err(N) = norm(y(tt) - sol(tt), inf);
    
end
fprintf('\n');

%% Plotting:

% Spy plot:
figure(6)
spy(A);
print -depsc2 example3_driscoll2_spy

% Error:
figure(3)
semilogy(NN, err(NN), '-', 'LineWidth', 3), shg
axis([0 NN(end) 1e-16, 1e1]); grid on
legend('New method', 'Driscoll (Chebfun)', 'Driscoll (integral)')
h = breakxaxis([60 980]);
set(h.leftAxes, 'XTick', [0 20 40 60]);
set(h.rightAxes, 'XTick', 1000, 'fontsize', 9.5);
legend(h.leftAxes, 'New method', 'Driscoll (Chebfun)', 'Driscoll (integral)')
print -depsc2 example3_err_all

alignfigs

