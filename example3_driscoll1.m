% Example 3 from VIDECT paper using Chebfun
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

% Chebfun
A = chebop(@(t,y) a^2*diff(y,2) + t*diff(y) + y + fred(@(s,t)K(t-s), y), d);
A.lbc = gam(1);
A.bc = @(t,y) integral(y) - gam(2);

err = [];
tt = linspace(0, T, 1000);
NN = [2:100 980:1000];
for N = [NN 60]
    fprintf('%d ', N);

    AN = matrix(A, N);
    
    rhs = [gam ; feval(f, chebpts(N,d,1))];
    
    y_ = AN\rhs;
    
    % Approximate infinity norm err
    y = chebfun(y_, d);
    err(N) = norm(y(tt) - sol(tt), inf);
    
end
fprintf('\n');

%% Plotting:

% Spy plot:
figure(5)
spy(A);
print -depsc2 example3_driscoll1_spy

% Error:
figure(3)
semilogy(NN, err(NN), '-', 'LineWidth', 3), shg
axis([0 NN(end) 1e-16, 1e1]); grid on

alignfigs