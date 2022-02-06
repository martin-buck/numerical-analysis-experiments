% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 2/5/22
%
% This script tests different fixed point iteration functions to find the
% root of a certain polynomial

f = @(x) x^2-3*x+2;
g1 = @(x) (x^2+2)/3;
g2 = @(x) sqrt(3*x-2);
g3 = @(x) 3-2/x;
g4 = @(x) (x^2-2)/(2*x-3);

g_list = {g1; g2; g3; g4};
titles = [" $g(x)=\frac{x^2+2}{3}, g'(2)=\frac{4}{3}$", " $g(x)=\sqrt(3x-2), g'(2)=\frac{3}{4}$", ...
    " $g(x)=3-\frac{2}{x}, g'(2)=\frac{1}{2}$", " $g(x)=\frac{x^2-2}{2x-3}, g'(2)=0$"];

% start the iteration close to the root of the polynomial x=2 since we are
% only testing convergence rates and divergence
x0 = 10;
root = 2;

% loop over the four functions and plot the iterates and the error
for i = 1:4
    figure
    g = g_list{i};
    iterates = zeros(1, 10);
    errors = zeros(1, 10);
    iterates(1) = g(x0);
    errors(1) = abs(iterates(1)-root);
    for j = 2:10
        % fixed point iteration for each of the four functions above
        iterates(j) = g(iterates(j-1));
        errors(j) = abs(iterates(j)-root);
        scatter(iterates(j-1), iterates(j), 'filled', 'r')
        hold on
    end
    
    x = linspace(0, 10);
    plot(x, x)
    grid on;
    xlabel('x', 'Interpreter', 'latex')
    ylabel('y', 'Interpreter', 'latex')
    title(strcat('Fixed Point Iterations, ', titles{i}), 'Interpreter', 'latex')
    hold off    
    
    figure
    scatter(1:10, errors, 'filled', 'b')
    grid on;
    xlabel('Iteration Number', 'Interpreter', 'latex');
    ylabel('Error ${e_k}$', 'Interpreter', 'latex')
    title(strcat('Error Ratio vs. Iteration Number, ', titles{i}), 'Interpreter', 'latex')
    hold on
    dx = .1; dy = .1;
    text(1:10+dx, errors+dy, cellstr(num2str(errors')))
end