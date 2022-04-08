% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 4/8/22
%
% This script implements the composite trapezoid rule for numerical
% integration and plots the error as a function of the grid spacing
% (regular grid) which should decrease like the square of the grid spacing

f1 = @(x) 1/(x+1);
f2 = @(x) exp(-x^2);
a = 0;
b = 1;
syms x;
approx_integrals = zeros(5, 1);
exact = int(sym(f2), x, a, b);

for i = 2:11
    N = 2^i;
    h = 1/N;
    approx_integrals(i-1) = trap(a, b, N, f2);
    scatter(h, approx_integrals(i-1), 'filled');
    hold on
end

set(gca, 'Xscale', 'log')
grid on
xlabel('Grid Spacing')
ylabel('Numerical Integral Approximation')
title('Composite Trapezoid Approximation vs. Mesh Spacing')

figure
for i = 2:11
    N = 2^i;
    h = 1/N;
    approx_integrals(i-1) = trap(a, b, N, f2);
    scatter(h, abs(approx_integrals(i-1)-exact), 'filled');
    hold on
end

set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
grid on
xlabel('Grid Spacing')
ylabel('Numerical Integral Approximation Error')
title('Composite Trapezoid Approximation Error vs. Mesh Spacing')

% this function takes an interval [a,b] as input and a parameter N
% corresponding to the number of panels in the composite trapezoid rule
function [I] = trap(a, b, N, f)
    % height of trapezoid
    h = (b-a)/N;
    I = 0;
    for i=1:N
        f1 = feval(f,(a+h*(i-1)));
        f2 = feval(f,(a+h*(i)));
        I = I + (h/2)*(f1 + f2);
    end
end
