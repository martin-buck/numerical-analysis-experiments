% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 3/9/22
%
% This script computes the solution to a system of non-linear equations
% using Broyden's method. This type of system is formed when deriving the
% two-pouint Gaussian Quadrature rule using the method of undetermined
% coefficients

% randomly initialize a starting point for the iteration
n = 4;
eps = 10e-4;
x0 = 2*rand(4, 1)-1 ;
% only one expensive Jacobian evaluation
B = gauss_quad_jacobian(x0);
f = @gauss_quad;
% do one step of Broyden so we can start the iteration based on making
% successive approximations small
s0 = -B\f(x0);
x1 = x0 + s0;
y1 = f(x1)-f(x0);
z = y1 - B*s0;
% rank-1 update to B. Should be a better approximation to the Jacobian!
B = B + (z*s0')/(s0'*s0);

sol_m = fsolve(f, x0);

sol_b = Broyden(x0, x1, B, f);

% Broyden's method
function [sol] = Broyden(x0, x1, B, f)
    xk = x0;
    xk1 = x1;
    
    while max(abs(xk1-xk)) > eps
        xk = xk1;
        sk = -B\f(xk);
        xk1 = xk + sk;
        yk = f(xk1)-f(xk);
        z = yk - B*sk;
        % rank-1 update to B. Should be a better approximation to the Jacobian!
        B = B + (z*sk')/(sk'*sk);
    end
    
    sol = xk1;
end

% The function whose zeroes we want to locate
function [f_out] = gauss_quad(x)
    f1 = x(3) + x(4) - 2;
    f2 = x(3)*x(1) + x(4)*x(2);
    f3 = x(3)*x(1)^2 + x(4)*x(2)^2 - 2/3;
    f4 = x(3)*x(1)^3 + x(4)*x(2)^3;
    f_out = [f1; f2; f3; f4];
end

% Evaluate the Jacobian. Only done once!
function [B] = gauss_quad_jacobian(x)
    B = [0 0 1 1; ...
        x(3) x(4) x(1) x(2); ...
        2*x(3)*x(1) 2*x(4)*x(2) x(1)^2 x(2)^2; ...
        3*x(3)*x(1)^2 3*x(4)*x(2)^2 x(1)^3 x(2)^3];
end