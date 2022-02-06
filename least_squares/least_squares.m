% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 2/4/22
%
% This script compares the solutions produced to the least squares problem
% by solving the normal equations (ill-conditioned) versus the SVD
% (well-conditioned)

m = 21;
n = 12;
t = linspace(0, 1, m);
x = ones(n, 1);

% evaluate a polynomial at 21 points in the unit interval, perturb it
% slightly, and then solve the least squares problem
[y, y_pert, A] = poly(t);
x_Chol = Chol_LS(A, y_pert);
x_SVD = SVD_LS(A, y_pert);

disp(max(x_Chol-x))
disp(max(x_SVD-x))

% create the Vandermonde matrix and the value of the polynomial at the
% different points in the interval and a slight perturbation of these
% values
function [y, y_pert, A] = poly(t)
    eps = 10^(-10); 
    m = 21;
    n = 12;
    A = zeros(m, n);
    x = ones(n, 1);
    
    for i=1:n
        A(:, i) = t.^(i-1);
    end
    y = A*x;
    y_pert = y + 2*eps*unifrnd(0, 1, m, 1);
end

% normal equations via the Cholesky factorization to produce a least
% squares solution
function [x] = Chol_LS(A, y_pert)
    z = A'*y_pert;
    R = chol(A'*A);
    c = R'\z;
    x = R\c;
end

% SVD decomposition to produce a least squares solution
function [x] = SVD_LS(A, y_pert)
    [U, S, V] = svd(A);
    y_U = U'*y_pert;
    x_V = S\y_U;
    x = V'\x_V;
end