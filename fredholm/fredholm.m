% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 4/27/22
%
% This script solves a Fredholm integral equation by using the composite
% trapezoidal rule to turn the integral into a linear system

f = @(x) ((x^2+1)^(3/2)-x^3)/3;
a = 0;
b = 1;

for i=3:2:9
    x = linspace(0, 1, i+1);
    f1 = arrayfun(f, x);
    A = fredholm_inner(i, a, b);
    u = A\f1';
    disp(cond(A))
    plot(x, u, '-o')
    hold on
end

plot(x, x, '-o')
grid on
xlabel('Interval [0,1]')
ylabel('Solution u(x) vs. y=x')
title('Numerical Solution to Fredholm Integral Equation')
legend('n=3', 'n=5', 'n=7', 'n=9', 'y=x')

for i=3:2:9
    figure
    legend_str = {};
    for j = 1:i
        x = linspace(0, 1, i+1);
        f1 = arrayfun(f, x);
        A = fredholm_inner(i, a, b);
        [U, S, V] = svd(A);
        u = zeros(i+1, 1);
        for k=1:j
            u = u + U(:, k)'*f1'*V(:, k)/S(k, k);
        end
        plot(x, u, '-o')
        hold on
        legend_str{j} = strcat('k=', num2str(k));
    end
    plot(x, x, '-o')
    grid on
    xlabel('Interval [0,1]')
    ylabel('Solution u(x) vs. y=x')
    title(strcat('SVD Numerical Solution to Fredholm Integral Equation, n=', num2str(k)))
    legend(legend_str);
end

function [A] = fredholm_inner(n, a, b)
    A = zeros(n+1, n+1);
    w = (b-a)/n;
    for i=1:n+1
        s = a + (i-1)/n;
        for j=1:n+1
            t = a + (j-1)/n;
            A(i,j) = w*(s^2+t^2)^(1/2);
        end
    end
    
    A(:, 1) = A(:, 1)/2;
    A(:, end) = A(:, end)/2;
end