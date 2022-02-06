% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 1/25/22
%
% This program plots the error of the finite difference and the centered
% difference approximation to the derivative for decreasing mesh spacings

% Iterate over mesh spacings and plot the error of the finite difference
for i = 0:16
    h = 10^(-i);
    error = abs(finite_diff(@tan, h, 1)-sec(1)^2);
    scatter(h, error, 'b', 'filled')
    hold on
end

set(gca,'xscale','log')
set(gca,'yscale','log')
grid on
title('Finite Difference Approximation to tan(x)')
xlabel('Grid Spacing')
ylabel('Error')

% Iterate over mesh spacings and plot the error of the central differences
figure
for i = 0:16
    h = 10^(-i);
    error = abs(centered_diff(@tan, h, 1)-sec(1)^2);
    scatter(h, error, 'b', 'filled')
    hold on
end

set(gca,'xscale','log')
set(gca,'yscale','log')
grid on
title('Centered Difference Approximation to tan(x)')
xlabel('Grid Spacing')
ylabel('Error')

% This code takes as input one of the 'Mathematical Functions' in MATLAB
% and returns the finite difference approximation to the derivative for the
% given value of x and mesh-spacing h
function [df] = finite_diff(f, h, x)
    df = (f(x+h)-f(x))/h;
end

% This code takes as input one of the 'Mathematical Functions' in MATLAB
% and returns the center difference approximation to the derivative for the
% given value of x and mesh-spacing h
function [df] = centered_diff(f, h, x)
    df = (f(x+h)-f(x-h))/(2*h);
end