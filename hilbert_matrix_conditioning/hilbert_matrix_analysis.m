% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 1/25/22
% 
% This script looks at how the condition number of the classical Hilbert
% matrix increases as the size of the matrix increases. The relationship
% between this condition number and the accuracy of the solution to a
% system of equations is explored

% Find the norm of the error and residual for each Hilbert matrix
for i = 2:16
    [H, b] = Hilbert(i);
    x = H\b;
    residual = b - H*x;
    max_norm_residual = max(residual);
    error_vec = ones(i, 1) - x;
    max_norm_error = max(error_vec);
       
    % Plot the norm of the error and residual and condition number
    scatter(i, log10(max_norm_residual), 'b', 'filled')
    hold on
    scatter(i, log10(max_norm_error), 'r', 'filled')
    hold on
    scatter(i, log10(cond(H)), 'k', 'filled')
    hold on
    legend('Residual', 'Error', 'K(H)', 'Location', 'northwest')
end

grid on
title('Residual, Error, and Condition Number for Hilbert System')
xlabel('Size of Hilbert Matrix')
ylabel('Log10 Infinity-norm of Error Vector and Condition Number')

% create the Hilbert matrix
function [H, b] = Hilbert(n)
    x = zeros(1, 2*n-1);
        for i = 1:2*n-1
            x(i) = 1/i;
        end
    H = zeros(n);
        for i = 1:n
            H(i, :) = x(i:i+n-1);
        end    
    b = H*ones(n, 1);
end